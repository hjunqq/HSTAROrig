#!/usr/bin/env python3
"""Small hand-written decks for reader branches that no real case reaches.

    extra_decks.py OUT_DIR

Each deck is cooks_membrane (2-D, 289 nodes, one group) with one switch set in 1.glb
and/or one auxiliary file written. All of them read to the end of global_data on the
current reader; use them with `sweep.py run --cases-root OUT_DIR` and trace_diff.py.
"""
import re
import shutil
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "cases/cases/cooks_membrane"
H = "text\ntext\n"  # the two title lines .nrt starts with

DECKS = {
    # (glb line, token) -> value ; file -> content
    "x_outind_oid": ({(6, 17): "-1"}, {"1.oid": "3\n 1 2 3\n"}),
    "x_nbackf_obsc": ({(6, 14): "1"}, {"1.obsc": (
        "header\n 1 1 2 3\n 5 1 2\n 5 6\n 0.5 0.5\n 7 2 1\n 7\n 1.0\n"
        " wstep 3\n 0.0 1.0 2.0\n 0.1 0.2 0.3\n 0.4 0.5 0.6\n")}),
    "x_nrt_99": ({}, {"1.nrt": H + " 1\n 2 99\n 10 2\n 11 12\n 0.5 0.5\n 20 1\n 21\n 1.0\n"}),
    "x_nrt_88_tral": ({}, {"1.nrt": H + " 1\n 2 88\n TRAL 10 2\n 11 12\n 0.5 0.5\n SKIP 20 1\n"}),
    "x_nrt_0_tral": ({}, {"1.nrt": H + " 1\n 2 0\n TRAL 10 2\n 11 12\n 0.5 0.5\n SKIP 20 1\n"}),
    "x_nrt_10_plate": ({}, {"1.nrt": H + " 1\n 1 10\n 1 10 11 1 2 1 2 0.5 0.5\n"}),
    "x_nrt_two_groups": ({}, {"1.nrt": H + " 2\n 1 99\n 10 1\n 11\n 1.0\n 1 0\n TRAL 30 2\n 31 32\n 0.3 0.7\n"}),
    "x_meshc1": ({(2, 11): "1"}, {"1.nrt": "text\n 1\n 10 2\n 11 12\n 0.5 0.5\ntext\ntext\n 1\n 1 0\n"}),
}


def main() -> int:
    out = Path(sys.argv[1])
    if out.exists():
        shutil.rmtree(out)
    for name, (edits, files) in DECKS.items():
        d = out / name
        shutil.copytree(SRC, d)
        lines = (d / "1.glb").read_bytes().split(b"\n")
        for (ln, tok), v in edits.items():
            code = lines[ln - 1]
            s, e = [m.span() for m in re.finditer(rb"[^\s,]+", code)][tok - 1]
            lines[ln - 1] = code[:s] + v.encode() + code[e:]
        (d / "1.glb").write_bytes(b"\n".join(lines))
        for f, text in files.items():
            (d / f).write_text(text)
    print(f"{len(DECKS)} decks in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
