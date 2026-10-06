#!/usr/bin/env python3
"""Make deck variants that push the reader into branches no real case reaches.

    mutate_glb.py OUT_ROOT --case NAME [...] --set RECORD:FIELD=V1,V2,...  [...]

RECORD names a .glb data record by the first title line above it that contains the
given word (case-insensitive), e.g. `npoin` or `type_problem`; FIELD is the 1-based
token on the data line right below that title. One variant per (case, record field,
value) is written to OUT_ROOT/<case>__<field>=<value>/ as a full copy of the case
with only that token changed.

The variants need not be physically meaningful or even readable: the oracle is that
two reader versions consume them identically, including failing on the same record.
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

CASES = Path(__file__).resolve().parents[2] / "cases/cases"


def find_data_line(lines: list[bytes], word: str) -> int:
    for i, ln in enumerate(lines[:-1]):
        if re.search(rb"(?i)\b" + re.escape(word.encode()) + rb"\b", ln):
            return i + 1
    raise LookupError(word)


def set_token(line: bytes, field: int, value: str) -> bytes:
    code, sep, comment = line.partition(b"!")
    spans = [m.span() for m in re.finditer(rb"[^\s,]+", code)]
    if field > len(spans):
        raise LookupError(f"field {field} beyond {len(spans)} tokens")
    s, e = spans[field - 1]
    return code[:s] + value.encode() + code[e:] + sep + comment


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--case", action="append", required=True)
    ap.add_argument("--set", action="append", required=True, help="WORD:FIELD=V1,V2")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    made = 0
    for case in a.case:
        glb = (CASES / case / "1.glb").read_bytes()
        lines = glb.split(b"\n")
        for spec in a.set:
            where, values = spec.split("=")
            word, field = where.split(":")
            try:
                k = find_data_line(lines, word)
                for v in values.split(","):
                    new = list(lines)
                    new[k] = set_token(lines[k], int(field), v)
                    dst = a.out / f"{case}__{word}{field}={v}"
                    if dst.exists():
                        shutil.rmtree(dst)
                    shutil.copytree(CASES / case, dst, symlinks=True)
                    (dst / "1.glb").write_bytes(b"\n".join(new))
                    made += 1
            except LookupError as exc:
                print(f"skip {case} {spec}: {exc}")
    print(f"{made} variants in {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
