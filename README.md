# Modernized BEM for Elastostatics (ELQUABE)

This repository contains a modernized Fortran 90 implementation of the **ELQUABE** program, originally developed by Brebbia and Dominguez in their classic textbook *Boundary Elements: An Introductory Course*.

The program uses the **Boundary Element Method (BEM)** with quadratic boundary elements to solve 2D (plane) elastostatic problems.

## 🚀 Modernization Features

The original Fortran 77 codebase has been fully refactored to meet modern high-performance computing standards:

1. **Free-Form & Modular Architecture**: The monolithic code was split into clean F90 modules (`data_mod`, `io_mod`, `assembly_mod`, `solver_mod`), eliminating all legacy `COMMON` blocks and enabling strict compiler type checking (`IMPLICIT NONE`).
2. **Dynamic Memory Allocation**: Static limitations (e.g., fixed dimensions `NX=100`, `NX1=150`) were completely removed. The program now reads the exact number of boundary elements and internal points from the input file and allocates arrays dynamically using `ALLOCATABLE`.
3. **OpenMP Parallelization**: The heavy matrix integration and assembly routines (`GHMATEQ` and `INTEREQ`) that build the $G$ and $H$ matrices and evaluate internal points are parallelized with OpenMP, mapping independent collocation and internal points to different CPU threads.
4. **LAPACK Solver Integration**: The native manual Gaussian elimination solver (`SLNPD`) was replaced with LAPACK's highly optimized `DGESV` (Double-precision Real General Linear Equation Solver) for superior performance and numerical stability in large systems.
5. **Double Precision Upgrade**: The solver inherently enforces **64-bit double precision** using compiler flags (`-fdefault-real-8`), promoting all `REAL` variables to `REAL*8` (`real64`), guaranteeing highly precise numerical integration.
6. **Python Pre/Post-Processing**: Python scripts (`preproc.py` and `posproc.py`) are provided to integrate with Gmsh, plot geometries, and validate numerical results against analytical elastostatic solutions.
7. **Colab Ready**: A Jupyter Notebook (`Run_elquabe_on_Colab.ipynb`) is included to easily run the entire pipeline in Google Colab.

## 📂 Project Structure

### Fortran Source Code
- `main.f90`: Main program entry point.
- `data_mod.f90`: Global variable declarations and dynamic array allocations.
- `assembly_mod.f90`: BEM influence matrices integration and assembly (`GHMATEQ`, `EXTINEQ`, `LOCINEQ`, `INTEREQ`, `SIGMAEQ`).
- `io_mod.f90`: Input/Output routines (`INPUTEQ`, `OUTPTEQ`).
- `solver_mod.f90`: System of equations solver wrapper for LAPACK (`DGESV`).

### Python Tools & Validation
- `preproc.py`: Reads Gmsh files (`.msh` via `meshio`), displays the geometry, and exports the simulation input data format (`eq12cb.dat`).
- `posproc.py`: Post-processing script that reads the exported boundary/internal displacements and stresses (`dpbou.dat`, `dpint.dat`, `stbou.dat`, `stint.dat`) and compares the numerical BEM solution against analytical plane elasticity solutions.
- `Run_elquabe_on_Colab.ipynb`: Google Colab notebook for running the preprocessing, compiling the BEM code, and executing the analytical validation entirely in the cloud.

### Build Scripts
- `Makefile` / `build.bat`: Build scripts for Windows/MSYS2 and Linux environments.

## ⚙️ Compilation & Build

### Requirements
- **GFortran** (GCC)
- **LAPACK & BLAS** libraries
- **OpenMP** (libgomp)

If you are using Windows with MSYS2, ensure the dependencies are installed:
```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-lapack mingw-w64-ucrt-x86_64-make
```

### Building the Executable
Navigate to the repository folder and simply run:
```bash
make
```
*(Or `mingw32-make` on MSYS2 Windows)*. This will invoke `gfortran` with `-fopenmp -O3 -fdefault-real-8` and link against `-llapack -lblas` to produce `elquabe_modern.exe`.

## 🏃 Usage & Validation Pipeline

1. **Preprocessing**: Create your geometry and mesh using Gmsh, and save it as `placa.msh` (or any `.msh`). Run the Python preprocessor to generate the standard input file:
   ```bash
   python preproc.py
   ```
   *Outputs: `eq12cb.dat` and `geometry.png`.*

2. **Execution**: Run the Fortran BEM solver:
   ```bash
   ./elquabe_modern.exe
   ```
   The solver will read `eq12cb.dat`, compute the elastostatic influence matrices across multiple threads, and solve the system via LAPACK (`DGESV`). It writes the formatted report to `eq12cb.out` and exports boundary/internal displacements and stresses to `dpbou.dat`, `dpint.dat`, `stbou.dat`, and `stint.dat`.

3. **Post-Processing & Validation**: Compare the BEM numerical results with the analytical solution for 2D plane elasticity:
   ```bash
   python posproc.py
   ```
   This script reads the numerical displacements and stresses, compares them against the theoretical elastostatic values, and exports the comparison plot to `result_plot.png`.

## 📖 References

* **Brebbia, C. A., & Dominguez, J.** (1992). *Boundary Elements: An Introductory Course*. Computational Mechanics Publications & McGraw-Hill.
