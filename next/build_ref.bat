@echo off
setlocal enabledelayedexpansion
REM Build the ORIGINAL (unmodified) HSTAR as reference for regression testing

set "HSTAR_DIR=G:\BB\hstarYLOrig\HSTAR"
set "BUILD_DIR=G:\BB\hstarYLOrig\next\build_ref"
set "LIB_DIR=G:\BB\hstarYLOrig\HSTAR\lib64"
set "MKL_INC=C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\include"
set "MKL_LIB=C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib"

call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat" intel64

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"

echo ============================================
echo  Building HSTAR Original (reference)
echo ============================================
echo.

cd /d "%BUILD_DIR%"

set FC=ifx
set FFLAGS=/nologo /Qmkl:parallel /I"%MKL_INC%"

echo [1/17] Vartype.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Vartype.f90" || goto :fail
echo [2/17] Array.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Array.f90" || goto :fail
echo [3/17] Elements.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Elements.f90" || goto :fail
echo [4/17] gidpost.F90
%FC% /c %FFLAGS% "%HSTAR_DIR%\gidpost.F90" || goto :fail
echo [5/17] Global.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Global.f90" || goto :fail
echo [6/17] Material.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Material.f90" || goto :fail
echo [7/17] meshfine.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\meshfine.f90" || goto :fail
echo [8/17] Prescrib.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Prescrib.f90" || goto :fail
echo [9/17] Load.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Load.f90" || goto :fail
echo [10/17] Output.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Output.f90" || goto :fail
echo [11/17] Solver.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Solver.f90" || goto :fail
echo [12/17] Stiff.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Stiff.f90" || goto :fail
echo [13/17] Residu.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Residu.f90" || goto :fail
echo [14/17] Temper.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Temper.f90" || goto :fail
echo [15/17] Level.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Level.f90" || goto :fail
echo [16/17] vsl_gauss_module.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\vsl_gauss_module.f90" || goto :fail
echo [17/17] Fem.f90 (original)
%FC% /c %FFLAGS% "%HSTAR_DIR%\Fem.f90" || goto :fail

echo.
echo Linking...
%FC% /nologo *.obj /Fe:hstar_ref.exe /link /NODEFAULTLIB:MSVCRT /LIBPATH:"%LIB_DIR%" /LIBPATH:"%MKL_LIB%" shlwapi.lib zlib.lib libhdf5_hl.lib libhdf5.lib gidpost.lib || goto :fail

echo.
echo BUILD SUCCEEDED: %BUILD_DIR%\hstar_ref.exe
goto :end

:fail
echo.
echo BUILD FAILED!
exit /b 1

:end
endlocal
