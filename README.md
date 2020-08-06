
# Geometric Volume-of-Fluid (GVOF) with Adaptive Mesh Refinement (AMR)

Contains test code to demonstrate standard Geometric volume-of-fluid techniques, such as the Piecewise linear interface construction (PLIC) technique
and its combination with AMR

The following dependencies must be available on your system:

- [AMReX](https://github.com/AMReX-Codes/amrex)


```
git clone https://github.com/gdeskos/AMR-GVOF.git
cd AMR-GVOF
mkdir build
cd build
cmake ..
make -j 4
```

## Usage

```
mpiexec -np <NPROCS> amr_gvof <input_file>
```
