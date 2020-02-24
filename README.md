
# Level-set and VOF for amr-wind

Contains test code to demonstrate multiphase capabilities of amr-wind

The following dependencies must be available on your system:

- [AMReX](https://github.com/AMReX-Codes/amrex)


```
git clone https://github.com/gdeskos/amr-offshore-utils.git
cd amr-offshore-utils
mkdir build
cd build
cp ../share/cmake_configure.sh 
# EDIT script as necessary
./cmake_configure.sh 
make 
```

## Usage

```
mpiexec -np <NPROCS> multpha-amr <input_file>
```
