# Parallel Spectral Deferred Correction

Benchmarck problem for the development of parallel Spectral Deferred Correction (SDC) method for time integration.
This provides scripts for the automated run and postprocessing of the the Galewsky test problem.

## Requirements

Uses base SWEET environnement that can be activated in the root directory :

```bash
$ source ./activate.sh
```
Need the following local installation in `local_software` :

```bash
$ cd local_software
$ ./install_miniconda.sh
$ ./install_scons.sh
$ ./install_numactl.sh
$ ./install_fftw3.sh
$ ./install_shtns.sh
$ ./install_shtns_python.sh
$ ./install_lapack.sh
```

### For GUI installation

Need the following local installation in `local_software` :

```bash
$ cd local_software
$ ./install_sdl2.sh
$ ./install_libfreetype.sh
```

And additionaly

```bash
$ sudo apt install pkg-config libgl-dev libxext-dev
```

### Using LLVM compiler

Installation on Ubuntu 22.04

```bash
# Clang compiler
sudo apt install clang-15 clangd-15 lldb-15 lld-15
# Additional c++ library (automatically included with gcc, not with clang)
sudo apt install libstdc++-12-dev
# Full openmp libraries
sudo apt install libomp-dev libomp5-15 libomp-15-dev
# And also gfortran separetly
sudo apt install libgfortran-12-dev
```

## Main scripts

- [1_create_jobs.py](./1_create_jobs.py) : create a running job and set all parameters using MULE. By default, it sets a parallel run using all available processors for space parallelization, but one can set the number of parallel processes like this (_e.g_ for 4 processes):
```bash
$ ./1_create_jobs.py 4
```
- [2_benchmark_compile.sh](./2_benchmark_compile.sh) : compile all the SWEET sources required for the test case run.
- [3_benchmark_run_all.sh](./3_benchmark_run_all.sh) : run all jobs instantiated using the `2_benchmark_compile.sh` script
- [4_postprocess.sh](./4_postprocess.sh) : postprocessing script (to be developped ...)