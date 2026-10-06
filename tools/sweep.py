#!/usr/bin/env python3
"""Run one HSTAR binary over the case library, each case in a private copy.

    sweep.py run --binary B --out DIR [--case NAME ...] [--timeout S] [--jobs N]
    sweep.py diff DIR_A DIR_B

`run` never touches cases/cases: every case directory is copied (symlinks kept as
links) into DIR/<case>, a missing `inp` gets the standard 4-line launcher, and the
binary runs single-threaded there. `report.json` records per case the status, exit
code, wall time, and sha256 of every file the run created or changed.

Status: OK (rc 0 and a finite 1.flavia.res with at least one block), NAN, NO_RESULT,
RC=<n>, TIMEOUT. Status is about the run, not about physical correctness.

`diff` compares two sweeps case by case: status, then 1.flavia.res values
(max |a-b|), then the hash of every other produced file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT.parent / "cases/cases"  # replaced by --cases-root
sys.path.insert(0, str(ROOT.parent / "hstar-evolution/tools"))
from yl_parse_flavia import parse  # noqa: E402

DEFAULT_INP = "input\n0 0 0 0 0 0  ! restart,relis,sysrelis,ADINA,Uopt_R,gamax\ninput\n1\n"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tree_state(d: Path) -> dict[str, tuple[int, int]]:
    return {str(p.relative_to(d)): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in d.rglob("*") if p.is_file() and not p.is_symlink()}


def flavia_status(res: Path) -> tuple[str, int]:
    if not res.exists() or res.stat().st_size == 0:
        return "NO_RESULT", 0
    try:
        parsed = parse(res)
    except Exception:  # noqa: BLE001 - an unparsable file is a failed run, not a crash here
        return "NO_RESULT", 0
    blocks = parsed["blocks"]
    if not blocks:
        return "NO_RESULT", 0
    for b in blocks:
        for row in b["rows"].values():
            if any(not math.isfinite(v) for v in row):
                return "NAN", len(blocks)
    return "OK", len(blocks)


def run_case(name: str, binary: Path, out: Path, timeout: int) -> dict:
    job = out / name
    shutil.copytree(CASES / name, job, symlinks=True)
    entry: dict = {"case": name}
    if not (job / "inp").exists():
        (job / "inp").write_text(DEFAULT_INP)
        entry["inp_generated"] = True
    before = tree_state(job)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", MKL_DYNAMIC="FALSE")
    t0 = time.monotonic()
    try:
        with open(job / "_stdout.log", "wb") as log:
            proc = subprocess.run([str(binary)], cwd=job, env=env, stdin=subprocess.DEVNULL,
                                  stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
        rc = proc.returncode
    except subprocess.TimeoutExpired:
        rc = None
    entry["seconds"] = round(time.monotonic() - t0, 2)
    entry["returncode"] = rc
    after = tree_state(job)
    entry["produced"] = {k: sha(job / k) for k, v in sorted(after.items())
                         if k != "_stdout.log" and before.get(k) != v}
    status, nblocks = flavia_status(job / "1.flavia.res")
    if "1.flavia.res" not in entry["produced"]:
        status, nblocks = "NO_RESULT", 0
    if rc is None:
        status = "TIMEOUT"
    elif rc != 0:
        status = f"RC={rc}"
    entry["status"], entry["result_blocks"] = status, nblocks
    tail = (job / "_stdout.log").read_bytes()[-2000:].decode("utf-8", "replace")
    entry["stdout_tail"] = tail.strip().splitlines()[-6:]
    return entry


def cmd_run(a) -> int:
    global CASES
    if a.cases_root:
        CASES = a.cases_root.resolve()
    binary = a.binary.resolve(strict=True)
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    names = a.case or sorted(p.name for p in CASES.iterdir()
                             if (p / "1.glb").exists() and not p.name.startswith("_"))
    report = {"binary": str(binary), "binary_sha256": sha(binary), "timeout": a.timeout,
              "cases": {}}
    with ThreadPoolExecutor(a.jobs) as pool:
        futs = {n: pool.submit(run_case, n, binary, out, a.timeout) for n in names}
        for n in names:
            e = futs[n].result()
            report["cases"][n] = e
            print(f"{e['status']:10s} {e['seconds']:8.1f}s  {n}", flush=True)
            (out / "report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False))
    counts: dict[str, int] = {}
    for e in report["cases"].values():
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    report["counts"] = counts
    (out / "report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False))
    print(counts)
    return 0


def max_diff(ra: Path, rb: Path) -> tuple[str, float]:
    try:
        pa, pb = parse(ra), parse(rb)
    except Exception as exc:  # noqa: BLE001
        return f"unparsable: {exc}", math.nan
    if [b["name"] for b in pa["blocks"]] != [b["name"] for b in pb["blocks"]]:
        return "block structure differs", math.nan
    m = 0.0
    for x, y in zip(pa["blocks"], pb["blocks"]):
        if x["rows"].keys() != y["rows"].keys():
            return f"node set differs in {x['name']}", math.nan
        for k, r in x["rows"].items():
            for u, v in zip(r, y["rows"][k]):
                if math.isfinite(u) and math.isfinite(v):
                    m = max(m, abs(u - v))
                elif not (math.isnan(u) and math.isnan(v)) and u != v:
                    return "non-finite mismatch", math.inf
    return ("identical" if m == 0.0 else "values differ"), m


CLOCK = re.compile(rb"\d\d:\d\d:\d\d")


def same_modulo_clock(pa: Path, pb: Path) -> bool:
    """The solver prints wall-clock stamps (`time: 01:59:49`) into 1.chk; ignore lines with one."""
    if not (pa.is_file() and pb.is_file()):
        return False
    def norm(p: Path) -> list[bytes]:
        return [ln for ln in p.read_bytes().splitlines() if not CLOCK.search(ln)]
    return norm(pa) == norm(pb)


def cmd_diff(a) -> int:
    A = json.loads((a.a / "report.json").read_text())["cases"]
    B = json.loads((a.b / "report.json").read_text())["cases"]
    same = 0
    for n in sorted(set(A) | set(B)):
        ea, eb = A.get(n), B.get(n)
        if ea is None or eb is None:
            print(f"{n}: only in {'A' if eb is None else 'B'}")
            continue
        notes = []
        if ea["status"] != eb["status"]:
            notes.append(f"status {ea['status']} -> {eb['status']}")
        fa, fb = ea["produced"], eb["produced"]
        if "1.flavia.res" in fa and "1.flavia.res" in fb and fa["1.flavia.res"] != fb["1.flavia.res"]:
            verdict, m = max_diff(a.a / n / "1.flavia.res", a.b / n / "1.flavia.res")
            if verdict != "identical":
                notes.append(f"flavia {verdict} max|d|={m:.3e}")
        other = sorted(k for k in set(fa) | set(fb)
                       if k != "1.flavia.res" and fa.get(k) != fb.get(k)
                       and not same_modulo_clock(a.a / n / k, a.b / n / k))
        if other:
            notes.append(f"{len(other)} other files differ: {', '.join(other[:6])}"
                         + (" ..." if len(other) > 6 else ""))
        if notes:
            print(f"{n} [{ea['status']}]: " + "; ".join(notes))
        else:
            same += 1
    print(f"identical cases: {same}/{len(set(A) & set(B))}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--binary", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--case", action="append")
    r.add_argument("--timeout", type=int, default=300)
    r.add_argument("--jobs", type=int, default=2)
    r.add_argument("--cases-root", type=Path, help="default: cases/cases")
    d = sub.add_parser("diff")
    d.add_argument("a", type=Path)
    d.add_argument("b", type=Path)
    a = ap.parse_args()
    return cmd_run(a) if a.cmd == "run" else cmd_diff(a)


if __name__ == "__main__":
    sys.exit(main())
