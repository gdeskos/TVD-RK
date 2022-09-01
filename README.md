
# TVD-RK with Adaptive Mesh Refinement (AMR)

Contains test code to demonstrate TOTAL VARIATION DIMINISHING RUNGE-KUTTA SCHEMES

The following dependencies must be available on your system:

- [AMReX](https://github.com/AMReX-Codes/amrex)


```
git clone https://github.com/gdeskos/AMR-GVOF.git
cd TVD-RK
mkdir build
cd build
cmake ..
make -j 4
```

## Usage

```
mpiexec -np <NPROCS> tvd_rk <input_file>
```
