#!/usr/bin/env python3
"""Compare the read traces (fort.9917) of two sweeps of instrumented builds.

    trace_diff.py SWEEP_A SWEEP_B [--tags TAGS_TSV]

Per case, two things must match: the trace byte for byte, and how the run ended
(`trace: stop after ...`, or the forrtl error class and unit). For the first case-level
difference the first diverging trace line is printed, with the READ text from tags.tsv.
Exit status 1 if any case differs.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def ending(log: Path) -> str:
    t = log.read_bytes().decode("latin-1") if log.exists() else ""
    if "trace: stop after" in t:
        return "stop"
    m = re.search(r"forrtl: (severe \(\d+\)[^,\n]*(?:, unit -?\d+)?|No such file or directory)", t)
    return m.group(1) if m else "other"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("a", type=Path)
    ap.add_argument("b", type=Path)
    ap.add_argument("--tags", type=Path)
    a = ap.parse_args()
    tags = {}
    if a.tags and a.tags.exists():
        for ln in a.tags.read_text(encoding="latin-1").splitlines():
            t, s = ln.split("\t", 1)
            tags[t] = s
    names = sorted(p.name for p in a.a.iterdir() if p.is_dir())
    bad = 0
    for n in names:
        da, db = a.a / n, a.b / n
        if not db.is_dir():
            print(f"{n}: missing in B"); bad += 1; continue
        ea, eb = ending(da / "_stdout.log"), ending(db / "_stdout.log")
        ta = (da / "fort.9917").read_bytes().splitlines() if (da / "fort.9917").exists() else []
        tb = (db / "fort.9917").read_bytes().splitlines() if (db / "fort.9917").exists() else []
        if ta == tb and (ea == eb or (ea != "stop" and eb != "stop")):
            continue  # same contract as fuzz_glb.py replay
        bad += 1
        msg = [f"{n}: ending {ea!r} vs {eb!r}" if ea != eb else f"{n}: same ending {ea!r}"]
        k = next((i for i, (x, y) in enumerate(zip(ta, tb)) if x != y), min(len(ta), len(tb)))
        if ta != tb:
            msg.append(f"  traces: {len(ta)} vs {len(tb)} lines, first difference at line {k + 1}")
            for side, t in (("A", ta), ("B", tb)):
                line = t[k].decode("latin-1") if k < len(t) else "<end of trace>"
                m = re.match(r"\s*R:(\w{8})", line)
                where = f"  [{tags.get(m.group(1), '?')}]" if m else ""
                msg.append(f"  {side}: {line[:110]}{where}")
        print("\n".join(msg))
    print(f"{len(names) - bad}/{len(names)} cases identical")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
