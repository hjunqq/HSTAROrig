#!/bin/bash
# Build hstarYLOrig for Linux using Intel ifx + MKL
set -e

SRCDIR="$(cd "$(dirname "$0")/HSTAR" && pwd)"
BUILDDIR="$(cd "$(dirname "$0")" && pwd)/build_linux"
MKL_INC="/opt/intel/oneapi/mkl/latest/include"
MKL_LIB="/opt/intel/oneapi/mkl/latest/lib/intel64"
FC=ifx

rm -rf "$BUILDDIR"
mkdir -p "$BUILDDIR"

# Compilation order (dependency order)
SRCS=(
    Vartype.f90
    Array.f90
    Elements.f90
    gidpost.F90
    vsl_gauss_module.f90
    Global.f90
    Material.f90
    meshfine.f90
    Load.f90
    Prescrib.f90
    Solver.f90
    Output.f90
    Temper.f90
    Stiff.f90
    Residu.f90
    Level.f90
    Fem.f90
)

echo "=== Building hstarYLOrig for Linux ==="
echo "Compiler: $FC"
echo ""

# Compile gidpost C stub
echo "Compiling gidpost C stub..."
STUB_SRC="$SRCDIR/../../_archive/src/gidpost_stub.c"
if [ ! -f "$STUB_SRC" ]; then
    echo "ERROR: gidpost_stub.c not found at $STUB_SRC"
    exit 1
fi
gcc -c "$STUB_SRC" -o "$BUILDDIR/gidpost_stub.o"
echo "  OK"

# Compile Fortran sources
OBJS=()
for f in "${SRCS[@]}"; do
    src="$SRCDIR/$f"
    obj="$BUILDDIR/${f%.*}.o"
    echo -n "  Compiling $f ..."
    $FC -c -O2 -module "$BUILDDIR" -I "$BUILDDIR" -I "$MKL_INC" \
        "$src" -o "$obj" 2>&1
    echo " OK"
    OBJS+=("$obj")
done

# Link
echo ""
echo "Linking..."
EXE="$SRCDIR/x64/Release/hstar"
$FC -O2 -qopenmp "${OBJS[@]}" "$BUILDDIR/gidpost_stub.o" -o "$EXE" \
    -L"$MKL_LIB" \
    -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
    -liomp5 -lpthread -lm -ldl

chmod +x "$EXE"
SIZE=$(du -h "$EXE" | cut -f1)
echo "  OK: $EXE ($SIZE)"
echo ""
echo "=== BUILD SUCCESSFUL ==="
