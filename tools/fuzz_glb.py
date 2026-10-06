#!/usr/bin/env python3
"""Grow a corpus of .glb variants that drive the input reader into many branches,
then replay it against another build and compare read traces.

    fuzz_glb.py grow   --binary TRACE_BIN --corpus DIR --case NAME [...] [--values 0,1,2,3,-1]
    fuzz_glb.py replay --binary TRACE_BIN --corpus DIR [--out DIR]
    fuzz_glb.py cover  --corpus DIR --tags TAGS_TSV

TRACE_BIN must be built from sources instrumented by trace_reads.py with
--stop-after, so a run only reads input and leaves fort.9917.

grow: for every integer token in a case's 1.glb and every value in --values, write the
variant and run it. When the run dies reading unit 1 (the .glb), a placeholder line of
ones is inserted where the failing READ would start, and the run is repeated, up to
--repairs times, so the variant gets past the branch it opened instead of masking
everything after it. The insertion point starts from "number of .glb READs that
succeeded" and tries a few lines around it (a READ may consume several lines); the
one that gets furthest is kept, and repair stops when none makes progress. Every run enters the corpus (see add()):
its 1.glb, its ending, and the hash of the trace it produced on TRACE_BIN; distinct
traces are stored once.

replay: run every corpus entry on another trace build and require the same trace
and the same ending. Exit status 1 on any difference.

cover: which instrumented READs the corpus reaches, and which it never does.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT.parent / "cases/cases"
PLACEHOLDER = b" ".join([b"1"] * 60)
TAG = re.compile(rb"^\s*R:(\w{8})")
ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", MKL_DYNAMIC="FALSE")


def ending(log: str) -> str:
    if "trace: stop after" in log:
        return "stop"
    m = re.search(r"forrtl: (severe \(\d+\)[^,\n]*(?:, unit -?\d+)?|No such file or directory)", log)
    if m:
        return m.group(1)
    return "timeout" if "<timeout>" in log else "other"


INPUT_EXT = {".glb", ".cor", ".ele", ".ftr", ".nrt", ".btl", ".obs", ".obsc", ".oid", ".resb"}


def _skip_big(d: str, names: list[str]) -> set[str]:
    out = set()
    for n in names:
        p = Path(d) / n
        if p.is_file() and not p.is_symlink() and p.suffix.lower() not in INPUT_EXT \
                and p.stat().st_size > 256 * 1024:
            out.add(n)
    return out


def run(binary: str, case: str, glb: bytes, timeout: int = 20) -> tuple[bytes, str]:
    with tempfile.TemporaryDirectory(dir=os.environ.get("FUZZ_TMP")) as d:
        job = Path(d) / "job"
        # Copying whole case directories dominates run time; large files that input
        # reading never touches (old results, restart dumps) are left out. If one of
        # them were read after all, base-vs-base replay would show it.
        shutil.copytree(CASES / case, job, symlinks=True, ignore=_skip_big)
        (job / "1.glb").write_bytes(glb)
        if not (job / "inp").exists():
            (job / "inp").write_text("input\n0 0 0 0 0 0\ninput\n1\n")
        try:
            p = subprocess.run([binary], cwd=job, env=ENV, stdin=subprocess.DEVNULL,
                               capture_output=True, timeout=timeout)
            log = (p.stdout + p.stderr).decode("latin-1")
        except subprocess.TimeoutExpired:
            log = "<timeout>"
        t = job / "fort.9917"
        return (t.read_bytes() if t.exists() else b""), ending(log)


def tag_seq(trace: bytes) -> list[str]:
    return [m.group(1).decode() for ln in trace.splitlines() if (m := TAG.match(ln))]


def int_tokens(glb: bytes) -> list[tuple[int, int, int, bytes]]:
    out = []
    for i, ln in enumerate(glb.split(b"\n")):
        code = ln.split(b"!")[0]
        for m in re.finditer(rb"[^\s,]+", code):
            if re.fullmatch(rb"-?\d+", m.group()):
                out.append((i, m.start(), m.end(), m.group()))
    return out


def grow_one(args) -> list[dict]:
    binary, case, values, repairs, glb_units = args
    base = (CASES / case / "1.glb").read_bytes()
    found = []
    for li, s, e, tok in int_tokens(base):
        for v in values:
            if v.encode() == tok:
                continue
            lines = base.split(b"\n")
            lines[li] = lines[li][:s] + v.encode() + lines[li][e:]
            glb = b"\n".join(lines)
            trace, end = run(binary, case, glb)
            found.append(entry(case, f"L{li + 1}:{tok.decode()}->{v}", 0, glb, trace, end))
            for r in range(1, repairs + 1):
                if not (end.startswith("severe") and end.endswith("unit 1")):
                    break
                # A READ whose list runs over several lines makes "READs so far" undercount
                # the lines consumed, so try a few insertion points and keep the one that
                # gets furthest (most READs; a clean stop wins).
                k = sum(1 for t in tag_seq(trace) if t in glb_units)
                best = None
                for off in (0, 1, 2, -1, 3):
                    lines = glb.split(b"\n")
                    if not 0 <= k + off <= len(lines):
                        continue
                    lines.insert(k + off, PLACEHOLDER)
                    cand = b"\n".join(lines)
                    t2, e2 = run(binary, case, cand)
                    score = (e2 == "stop", len(tag_seq(t2)))
                    if best is None or score > best[0]:
                        best = (score, cand, t2, e2)
                if best is None or best[0][1] <= len(tag_seq(trace)) and not best[0][0]:
                    break  # no insertion point made progress
                _, glb, trace, end = best
                found.append(entry(case, f"L{li + 1}:{tok.decode()}->{v}", r, glb, trace, end))
    return found


def entry(case, mutation, repairs, glb, trace, end) -> dict:
    seq = tag_seq(trace)
    key = hashlib.sha1((",".join(seq) + "|" + end).encode()).hexdigest()[:16]
    return {"case": case, "mutation": mutation, "repairs": repairs, "ending": end,
            "key": key, "glb": glb, "trace": trace, "tags": sorted(set(seq))}


def load_units(tags_tsv: Path) -> set[str]:
    out = set()
    for ln in tags_tsv.read_text(encoding="latin-1").splitlines():
        t, s = ln.split("\t", 1)
        if re.match(r"(?i)read\s*\(\s*gunit\b", s):
            out.add(t)
    return out


def cmd_grow(a) -> int:
    corpus = a.corpus
    corpus.mkdir(parents=True, exist_ok=True)
    manifest_p = corpus / "manifest.json"
    manifest = json.loads(manifest_p.read_text()) if manifest_p.exists() else {}
    units = load_units(a.tags)
    jobs = [(str(a.binary.resolve()), c, a.values.split(","), a.repairs, units) for c in a.case]
    # the real cases themselves are corpus entries too
    for c in a.case:
        glb = (CASES / c / "1.glb").read_bytes()
        trace, end = run(str(a.binary.resolve()), c, glb)
        add(corpus, manifest, entry(c, "original", 0, glb, trace, end))
    with ProcessPoolExecutor(a.jobs) as pool:
        for c, found in zip(a.case, pool.map(grow_one, jobs)):
            before = len(manifest)
            for e in found:
                add(corpus, manifest, e)
            ends = {}
            for e in found:
                ends[e["ending"]] = ends.get(e["ending"], 0) + 1
            print(f"{c}: {len(found)} runs, {len(manifest) - before} new entries, endings {ends}",
                  flush=True)
            manifest_p.write_text(json.dumps(manifest, indent=1))
    manifest_p.write_text(json.dumps(manifest, indent=1))
    print(f"corpus: {len(manifest)} entries")
    return 0


def add(corpus: Path, manifest: dict, e: dict) -> None:
    """Every run is kept: deduplicating by the reference trace would throw away exactly
    the inputs that tell two reader versions apart (nlayer=3 traces like nlayer=0 under
    `nlayer==2`, but not under `nlayer>=2`). Decks are small; full traces are stored once
    per distinct trace, entries refer to them by hash."""
    glb_sha = hashlib.sha1(e["glb"]).hexdigest()[:16]
    eid = f"{e['case']}:{glb_sha}"
    if eid in manifest:
        return
    tsha = hashlib.sha1(e["trace"]).hexdigest()[:16]
    (corpus / "decks").mkdir(exist_ok=True)
    (corpus / "traces").mkdir(exist_ok=True)
    (corpus / "decks" / f"{glb_sha}.glb").write_bytes(e["glb"])
    tp = corpus / "traces" / tsha
    if not tp.exists():
        tp.write_bytes(e["trace"])
    manifest[eid] = {"case": e["case"], "deck": glb_sha, "trace": tsha,
                     **{k: e[k] for k in ("mutation", "repairs", "ending", "tags")}}


def replay_one(args):
    binary, corpus, key, case, deck = args
    glb = (Path(corpus) / "decks" / f"{deck}.glb").read_bytes()
    trace, end = run(binary, case, glb)
    return key, trace, end


def cmd_replay(a) -> int:
    manifest = json.loads((a.corpus / "manifest.json").read_text())
    if a.case:
        manifest = {k: m for k, m in manifest.items() if m["case"] in a.case}
    jobs = [(str(a.binary.resolve()), str(a.corpus), k, m["case"], m["deck"])
            for k, m in manifest.items()]
    bad = 0
    seen: dict[str, int] = {}
    with ProcessPoolExecutor(a.jobs) as pool:
        for key, trace, end in pool.map(replay_one, jobs, chunksize=8):
            m = manifest[key]
            if a.collect:
                for t in set(tag_seq(trace)):
                    seen[t] = seen.get(t, 0) + 1
                continue
            ref = (a.corpus / "traces" / m["trace"]).read_bytes()
            if trace != ref or end != m["ending"]:
                bad += 1
                if bad <= 10:
                    ra, rb = ref.splitlines(), trace.splitlines()
                    k = next((i for i, (x, y) in enumerate(zip(ra, rb)) if x != y), min(len(ra), len(rb)))
                    print(f"{key} {m['case']} {m['mutation']} (+{m['repairs']} repairs): "
                          f"ending {m['ending']!r} -> {end!r}; trace {len(ra)} vs {len(rb)} lines, "
                          f"first difference at line {k + 1}")
                    for side, t in (("ref", ra), ("new", rb)):
                        print(f"   {side}: {t[k][:100].decode('latin-1') if k < len(t) else '<end>'}")
    if a.collect:
        a.collect.write_text(json.dumps(seen, indent=1))
        print(f"{len(seen)} READ tags reached; written to {a.collect}")
        return 0
    print(f"{len(manifest) - bad}/{len(manifest)} corpus entries identical")
    return 1 if bad else 0


def cmd_cover(a) -> int:
    manifest = json.loads((a.corpus / "manifest.json").read_text())
    tags = dict(ln.split("\t", 1) for ln in a.tags.read_text(encoding="latin-1").splitlines())
    real = set().union(*(set(m["tags"]) for m in manifest.values() if m["mutation"] == "original"))
    allc = set().union(*(set(m["tags"]) for m in manifest.values()))
    print(f"distinct READ statements instrumented: {len(tags)}")
    print(f"reached by the original cases: {len(real & tags.keys())}")
    print(f"reached by the corpus:         {len(allc & tags.keys())}")
    print("never reached:")
    for t in sorted(tags.keys() - allc, key=lambda t: tags[t]):
        print(f"  {t}  {tags[t]}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("grow")
    g.add_argument("--binary", type=Path, required=True)
    g.add_argument("--corpus", type=Path, required=True)
    g.add_argument("--tags", type=Path, required=True)
    g.add_argument("--case", action="append", required=True)
    g.add_argument("--values", default="0,1,2,3,-1")
    g.add_argument("--repairs", type=int, default=6)
    g.add_argument("--jobs", type=int, default=2)
    r = sub.add_parser("replay")
    r.add_argument("--binary", type=Path, required=True)
    r.add_argument("--corpus", type=Path, required=True)
    r.add_argument("--jobs", type=int, default=2)
    r.add_argument("--case", action="append", help="only entries derived from these cases")
    r.add_argument("--collect", type=Path, help="only record which READ tags run, into this file")
    c = sub.add_parser("cover")
    c.add_argument("--corpus", type=Path, required=True)
    c.add_argument("--tags", type=Path, required=True)
    a = ap.parse_args()
    return {"grow": cmd_grow, "replay": cmd_replay, "cover": cmd_cover}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
