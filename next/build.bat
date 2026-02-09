@echo off
setlocal enabledelayedexpansion
REM Build script for HSTAR Phase 1 (decomposed Fem.f90)

REM Save paths before calling env scripts (they can mess up %~dp0)
set "HSTAR_DIR=G:\BB\hstarYLOrig\HSTAR"
set "NEXT_DIR=G:\BB\hstarYLOrig\next"
set "BUILD_DIR=G:\BB\hstarYLOrig\next\build"
set "LIB_DIR=G:\BB\hstarYLOrig\HSTAR\lib64"

REM Set up environment
call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\env\vars.bat" intel64

REM Create build directory
if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"

echo ============================================
echo  Building HSTAR Phase 1 (decomposed)
echo ============================================
echo.

cd /d "%BUILD_DIR%"

set FC=ifx
set "MKL_INC=C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\include"
set "MKL_LIB=C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib"
set FFLAGS=/nologo /debug:full /Od /warn:interfaces /traceback /check:bounds /Qmkl:parallel /I"%NEXT_DIR%" /I"%MKL_INC%"

REM Compile each file in dependency order
echo [1/18] Vartype.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Vartype.f90" || goto :fail
echo [2/18] Array.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Array.f90" || goto :fail
echo [3/18] Elements.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Elements.f90" || goto :fail
echo [4/18] gidpost.F90
%FC% /c %FFLAGS% "%HSTAR_DIR%\gidpost.F90" || goto :fail
echo [5/18] Global.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Global.f90" || goto :fail
echo [6/18] Material.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Material.f90" || goto :fail
echo [7/18] meshfine.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\meshfine.f90" || goto :fail
echo [8/18] Prescrib.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Prescrib.f90" || goto :fail
echo [9/18] Load.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Load.f90" || goto :fail
echo [10/18] Output.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Output.f90" || goto :fail
echo [11/18] Solver.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Solver.f90" || goto :fail
echo [12/18] Stiff.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Stiff.f90" || goto :fail
echo [13/18] Residu.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Residu.f90" || goto :fail
echo [14/18] Temper.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Temper.f90" || goto :fail
echo [15/18] Level.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\Level.f90" || goto :fail
echo [16/18] vsl_gauss_module.f90
%FC% /c %FFLAGS% "%HSTAR_DIR%\vsl_gauss_module.f90" || goto :fail
echo [17/18] fem_module.f90
%FC% /c %FFLAGS% "%NEXT_DIR%\fem_module.f90" || goto :fail
echo [18/18] Fem.f90
%FC% /c %FFLAGS% "%NEXT_DIR%\Fem.f90" || goto :fail

echo.
echo Linking...
%FC% /nologo *.obj /Fe:hstar_phase1.exe /link /NODEFAULTLIB:MSVCRT /LIBPATH:"%LIB_DIR%" /LIBPATH:"%MKL_LIB%" shlwapi.lib zlib.lib libhdf5_hl.lib libhdf5.lib gidpost.lib || goto :fail

echo.
echo BUILD SUCCEEDED: %BUILD_DIR%\hstar_phase1.exe
echo.
goto :end

:fail
echo.
echo BUILD FAILED!
exit /b 1

:end
endlocal
