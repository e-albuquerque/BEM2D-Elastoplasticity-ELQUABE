@echo off
REM ============================================================================
REM Build script for Modernized ELQUABE on Windows (MSYS2 / MinGW-w64 UCRT64)
REM ============================================================================

if exist "C:\msys64\ucrt64\bin" (
    set "PATH=C:\msys64\ucrt64\bin;%PATH%"
)

echo Compiling Modernized ELQUABE...
gfortran -fopenmp -O3 -fdefault-real-8 -Wall -Wextra ^
    data_mod.f90 io_mod.f90 assembly_mod.f90 solver_mod.f90 main.f90 ^
    -o elquabe_modern.exe -llapack -lblas

if %ERRORLEVEL% EQU 0 (
    echo Build successful: elquabe_modern.exe
) else (
    echo Build failed with error code %ERRORLEVEL%.
    exit /b %ERRORLEVEL%
)
