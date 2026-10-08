#!/bin/bash
# Debug build: -O0 -g -traceback -check bounds,pointers → x64/Release/hstar_debug
# Never touches the installed release binary.
set -e

SRCDIR="$(cd "$(dirname "$0")/HSTAR" && pwd)"
BUILDDIR="$(cd "$(dirname "$0")" && pwd)/build_linux_debug"
MKL_INC="/opt/intel/oneapi/mkl/latest/include"
MKL_LIB="/opt/intel/oneapi/mkl/latest/lib/intel64"
FC=ifx
FFLAGS="-O0 -g -traceback -check bounds,pointers"

rm -rf "$BUILDDIR"
mkdir -p "$BUILDDIR"

SRCS=(Vartype.f90 Array.f90 Elements.f90 gidpost.F90 vsl_gauss_module.f90
      Global.f90 Material.f90 meshfine.f90 Load.f90 Prescrib.f90 Solver.f90
      Output.f90 Temper.f90 Stiff.f90 Residu.f90 Level.f90 Fem.f90)

gcc -c "$SRCDIR/../../_archive/src/gidpost_stub.c" -o "$BUILDDIR/gidpost_stub.o"

OBJS=()
for f in "${SRCS[@]}"; do
    obj="$BUILDDIR/${f%.*}.o"
    echo -n "  $f ..."
    $FC -c $FFLAGS -module "$BUILDDIR" -I "$BUILDDIR" -I "$MKL_INC" \
        "$SRCDIR/$f" -o "$obj"
    echo " OK"
    OBJS+=("$obj")
done

EXE="$SRCDIR/x64/Release/hstar_debug"
$FC $FFLAGS -qopenmp "${OBJS[@]}" "$BUILDDIR/gidpost_stub.o" -o "$EXE" \
    -L"$MKL_LIB" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
    -liomp5 -lpthread -lm -ldl
chmod +x "$EXE"
echo "OK: $EXE"
