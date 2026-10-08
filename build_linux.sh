#!/bin/bash
# Build HSTAR for Linux with Intel ifx + MKL.
#   ./build_linux.sh [release|debug] [OUTDIR]
# Requires the oneAPI environment (MKLROOT set), e.g.
#   source /home/hatch/intel/oneapi/setvars.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SRCDIR="${HSTAR_SRC:-$ROOT/HSTAR}"   # override to build an instrumented copy
PROFILE="${1:-release}"
BUILDDIR="${2:-$ROOT/build/$PROFILE}"
STUB_SRC="${GIDPOST_STUB:-$ROOT/../_archive/src/gidpost_stub.c}"
: "${MKLROOT:?MKLROOT not set - source oneAPI setvars.sh first}"

case "$PROFILE" in
    release) FFLAGS=(-O2) ;;
    trace)   FFLAGS=(-O0) ;;   # read-trace builds: fast to compile, values identical
    debug)   FFLAGS=(-O0 -g -traceback -check bounds,pointers,uninit -fpe0) ;;
    *) echo "unknown profile: $PROFILE" >&2; exit 2 ;;
esac

# Dependency order
SRCS=(Vartype.f90 Array.f90 Elements.f90 gidpost.F90 vsl_gauss_module.f90
      yl_authoring_toml.f90 NativeInput.f90 Global.f90 Material.f90 meshfine.f90 Load.f90 Prescrib.f90 Solver.f90
      Output.f90 Temper.f90 Stiff.f90 Residu.f90 Level.f90 Fem.f90)

rm -rf "$BUILDDIR"; mkdir -p "$BUILDDIR"
[ -f "$STUB_SRC" ] || { echo "gidpost stub not found: $STUB_SRC" >&2; exit 1; }
icx -c -O2 "$STUB_SRC" -o "$BUILDDIR/gidpost_stub.o"

OBJS=()
for f in "${SRCS[@]}"; do
    obj="$BUILDDIR/${f%.*}.o"
    echo "  $f"
    ifx -c "${FFLAGS[@]}" -qopenmp -module "$BUILDDIR" -I "$BUILDDIR" \
        -I "$SRCDIR" -I "$MKLROOT/include" "$SRCDIR/$f" -o "$obj"
    OBJS+=("$obj")
done

ifx "${FFLAGS[@]}" -qopenmp "${OBJS[@]}" "$BUILDDIR/gidpost_stub.o" \
    -o "$BUILDDIR/hstar" -qmkl=parallel -lpthread -lm -ldl
echo "built: $BUILDDIR/hstar"
