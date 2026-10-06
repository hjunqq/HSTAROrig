#!/bin/bash
# Prove a freshly built solver reproduces the installed one, before replacing it.
#
# Rebuilding the authoritative solver is only safe if the new binary gives the
# same answer on a case we already trust. This runs one known-good deck
# (train01_gravdam_static, part of the 14/15 suite) through BOTH binaries in
# clean temp dirs and compares 1.flavia.res numerically.
#
#   ./verify_build.sh <candidate-binary> [case-dir]
#
# Exit 0 = the candidate matches. Anything else = do not install.
set -e

NEW="${1:?usage: verify_build.sh <candidate-binary> [case-dir]}"
CASE="${2:-/home/huijun/HSTAR_Next/cases/cases/train01_gravdam_static}"
CUR="$(cd "$(dirname "$0")" && pwd)/x64/Release/hstar"
TOL="${VERIFY_TOL:-1e-9}"

[ -x "$NEW" ] || { echo "candidate not executable: $NEW"; exit 2; }
[ -d "$CASE" ] || { echo "case dir not found: $CASE"; exit 2; }

WORK=$(mktemp -d /tmp/hstar-verify-XXXX)
trap 'rm -rf "$WORK"' EXIT

for side in cur new; do
  mkdir -p "$WORK/$side"
  # Inputs only — a stale result would let the comparison pass without solving.
  for f in "$CASE"/*; do
    b=$(basename "$f")
    case "$b" in
      1.flavia.res|1.flavia.msh|1.chk|1.res|1.gpv|1.act|1.ftf|1.dis|run.log|*.vtp|*.vtu|*.pvd|result_b*) ;;
      *) cp -a "$f" "$WORK/$side/" 2>/dev/null || true ;;
    esac
  done
done

run_one() {  # $1 = binary, $2 = dir
  ( cd "$2" && HSTAR_EXE="$1" timeout 1800 python3 \
      /home/huijun/HSTAR_Next/fem-chat/harness/_solve_deck.py . >solve.log 2>&1 ) || true
}

echo "  case: $(basename "$CASE")"
echo "  installed: $CUR"
run_one "$CUR" "$WORK/cur"
echo "  candidate: $NEW"
run_one "$NEW" "$WORK/new"

python3 - "$WORK/cur/1.flavia.res" "$WORK/new/1.flavia.res" "$TOL" <<'PY'
import sys
a, b, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])

def nums(p):
    out = []
    try:
        for line in open(p, errors="replace"):
            f = line.split()
            if len(f) < 2:
                continue
            try:
                int(f[0]); out.extend(float(x) for x in f[1:])
            except ValueError:
                continue
    except OSError:
        return None
    return out

x, y = nums(a), nums(b)
if not x or not y:
    print(f"  FAIL: missing or empty result (installed={x and len(x)}, candidate={y and len(y)})")
    sys.exit(1)
if len(x) != len(y):
    print(f"  FAIL: value count differs — installed {len(x)}, candidate {len(y)}")
    sys.exit(1)
scale = max(abs(v) for v in x) or 1.0
worst = max(abs(p - q) for p, q in zip(x, y)) / scale
print(f"  {len(x)} values, worst relative difference {worst:.3E} (tol {tol:.0E})")
sys.exit(0 if worst <= tol else 1)
PY
