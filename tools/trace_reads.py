#!/usr/bin/env python3
"""Make an instrumented copy of the HSTAR sources that traces every READ in chosen routines.

    trace_reads.py SRC_DIR OUT_DIR --routine Global.f90:global_data [--routine ...]
                   [--stop-after global_data]

For every READ statement inside a selected routine, the copy gets, right after the
read, a list-directed WRITE to unit 9917 (file `fort.9917` in the run directory):

    'R:<tag>', <the read's own input list>

<tag> is 8 hex digits of sha1 over the read statement's normalised text, so the same
statement gets the same tag wherever a refactor moves it. `tags.tsv` in OUT_DIR maps
tags back to text. One-line `if (c) read ...` becomes an IF block so the trace sits
under the same condition.

`--stop-after R` inserts `stop 'trace: stop after R'` right after `call R` in Fem.f90,
so a run ends when input reading for R is done: reading-only, no solve.

The oracle this serves: two source versions, instrumented by this same script, must
write identical fort.9917 for the same input - same records consumed in the same
order, same values assigned, and a failure on the same record. It relies on one
refactor rule: READ statement text is never edited, only moved.

Files are handled as bytes (latin-1 round trip): several HSTAR sources are GBK.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
from pathlib import Path

UNIT = 9917
SUB_START = re.compile(r"^\s*(?:recursive\s+)?subroutine\s+(\w+)", re.I)
SUB_END = re.compile(r"^\s*end\s*subroutine\b", re.I)


def strip_comment(line: str) -> str:
    """Code part of a free-form line (drop `!` comments outside string literals)."""
    q = None
    for i, ch in enumerate(line):
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "!":
            return line[:i]
    return line


def match_paren(s: str, i: int) -> int:
    """Index just past the ')' matching the '(' at s[i]."""
    depth, q = 0, None
    for j in range(i, len(s)):
        ch = s[j]
        if q:
            if ch == q:
                q = None
            continue
        if ch in "'\"":
            q = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return j + 1
    raise ValueError(f"unbalanced parentheses: {s}")


def wrap(stmt: str, indent: str) -> list[str]:
    """Emit one statement as free-form lines of <= ~100 chars, breaking after commas."""
    out, cur, q, depth = [], "", None, 0
    for ch in stmt:
        cur += ch
        if q:
            if ch == q:
                q = None
            continue
        if ch in "'\"":
            q = ch
        elif ch == "," and len(cur) > 90:
            out.append(cur + "&")
            cur = "    &"
    out.append(cur)
    return [indent + x for x in out]


def parse_read(code: str):
    """Split a logical statement into (label, if_cond, read_stmt, input_list) or None."""
    m = re.match(r"\s*(\d+\s+)?", code)
    label = (m.group(1) or "").strip()
    rest = code[m.end():]
    cond = None
    if re.match(r"if\s*\(", rest, re.I):
        p = rest.index("(")
        e = match_paren(rest, p)
        tail = rest[e:].lstrip()
        if not re.match(r"read\s*\(", tail, re.I):
            return None
        cond = rest[p:e]
        rest = tail
    if not re.match(r"read\s*\(", rest, re.I):
        return None
    p = rest.index("(")
    e = match_paren(rest, p)
    return label, cond, rest.strip(), rest[e:].strip()


def instrument(path: Path, routines: set[str], tags: dict[str, str], site: bool = False) -> int:
    lines = path.read_bytes().decode("latin-1").split("\n")
    out, i, current, count = [], 0, None, 0
    while i < len(lines):
        line = lines[i]
        m = SUB_START.match(strip_comment(line))
        if m:
            current = m.group(1).lower()
        elif SUB_END.match(strip_comment(line)):
            current = None
        if current not in routines:
            out.append(line)
            i += 1
            continue
        # gather one logical statement (continuations end in '&')
        j, parts = i, []
        while True:
            code = strip_comment(lines[j]).rstrip()
            if code.endswith("&"):
                parts.append(code[:-1])
                j += 1
                continue
            parts.append(code)
            break
        joined = "".join(re.sub(r"^\s*&", "", p) if k else p for k, p in enumerate(parts))
        parsed = parse_read(joined) if ";" not in joined else None
        if not parsed:
            out.append(line)
            i += 1
            continue
        label, cond, read_stmt, ilist = parsed
        norm = re.sub(r"\s+", "", read_stmt).lower()
        if site:  # coverage builds: one tag per source site, not per statement text
            norm += f"#{path.name}:{i + 1}"
        tag = hashlib.sha1(norm.encode("latin-1")).hexdigest()[:8]
        tags[tag] = read_stmt + (f"\t{path.name}:{i + 1}" if site else "")
        indent = re.match(r"\s*", lines[i]).group(0) or "    "
        trace = f"write({UNIT},*) 'R:{tag}'" + (f", {ilist}" if ilist else "")
        if cond is None:
            out.extend(lines[i:j + 1])
            out.extend(wrap(trace, indent))
        else:
            head = (label + " " if label else "") + f"if {cond} then"
            out.extend(wrap(head, indent))
            out.extend(wrap(read_stmt, indent + "    "))
            out.extend(wrap(trace, indent + "    "))
            out.append(indent + "endif")
        count += 1
        i = j + 1
    path.write_bytes("\n".join(out).encode("latin-1"))
    return count


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--routine", action="append", required=True, help="FILE:routine")
    ap.add_argument("--stop-after")
    ap.add_argument("--site", action="store_true",
                    help="tag each READ by its source line (coverage only; not stable under refactoring)")
    a = ap.parse_args()
    if a.out.exists():
        shutil.rmtree(a.out)
    shutil.copytree(a.src, a.out, ignore=shutil.ignore_patterns("x64", "*.bak*", "lib64"))
    by_file: dict[str, set[str]] = {}
    for r in a.routine:
        f, name = r.split(":")
        by_file.setdefault(f, set()).add(name.lower())
    tags: dict[str, str] = {}
    for f, names in by_file.items():
        n = instrument(a.out / f, names, tags, a.site)
        print(f"{f}: {n} reads instrumented in {sorted(names)}")
    if a.stop_after:
        fem = a.out / "Fem.f90"
        text = fem.read_bytes().decode("latin-1")
        pat = re.compile(rf"^(\s*)call\s+{a.stop_after}\b.*$", re.I | re.M)
        if len(pat.findall(text)) != 1:
            raise SystemExit(f"expected exactly one `call {a.stop_after}` in Fem.f90")
        text = pat.sub(lambda m: m.group(0) + f"\n{m.group(1)}stop 'trace: stop after {a.stop_after}'", text)
        fem.write_bytes(text.encode("latin-1"))
    (a.out / "tags.tsv").write_text("".join(f"{t}\t{s}\n" for t, s in sorted(tags.items())),
                                    encoding="latin-1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
