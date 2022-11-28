Before trying out anything, please have a detailed look into the INSTALL.md instructions about how to cope with the required libraries!

# SWEET

Shallow Water Equation Environment for Tests, Awesome!

Please checkout the website for more information
https://schreiberx.github.io/sweetsite

This library supports various kinds of simulations (not only Shallow Water)
on 2D surfaces. The currently supported surfaces are

 * the 2D torus (bi-periodic plane) and
 * the sphere

For discretization, we use the double FFT for the plane and Spherical Harmonics
on the sphere.


## Community

There's a discord channel for the SWEET community.
Drop an email to `martin.schreiber@univ-grenoble-alpes.fr` to get an invite.


## File structure

 * `archive/`: Old stuff which should be still visible when checking out repository
 * `benchmarks_sphere/`:	Scripts for running benchmarks (numerical and HPC) on the sphere
 * `benchmarks_generic/`:	Benchmarks without plane or sphere spatial discretization
 * `benchmarks_plane/`:	Scripts for running benchmarks (numerical and HPC) on the plane
 * `build/`:		Build directory, note, that also /tmp is used for storing object files
 * `data/`:		Data such as OpenGL shaders
 * `doc/`:		Documentation
 * `local_software/`:	3rd party software with installer scripts.
			If you get compile errors, please checkout this directory
 * `mule/`:		Compile, run and parallelization (PRC) framework
 * `mule_local/`:	Software-specific PRC scripts
 * `src/`:		Source folder
 * `src/programs`:	Example programs
 * `src/unit_tests`:	Unit tests to validate SWEET
 * `src/include/libgl`:	Some visualization helper tools
 * `src/include/libmath`:	Some mathematical routines
 * `src/include/sweet`:	DUDE! That's SWEET!
 * `src/include/...`:		Plenty of other stuff

 * `run_tests_compile_all.sh`:	Test script to check compilation of programs
 * `run_tests_validation.sh`:		Test script to validate computations
 * `run_valgrind.sh`:		Convenient script to call valgrind with preset options (Useful for debugging)
 * `SConstruct`:	Required for compilation (scons)



## SWEET Environment

First of all, setup the environment variables correctly (see also INSTALL file).


```
$ source ./activate.sh
```

## Compilation on personal workstation

After loading the environment, the compilation is done by calling 'scons' which is a
makefile alternative and offers the flexibility of Python.

```
$ scons --program=swe_plane
```

The possible compilation (maybe not complete) options are visible via

```
$ scons --help
```

The program to be compiled can be specified via -program=...
with the list of programs given in src/examples



## First steps

As a starter, you execute the benchmark in 
`benchmarks_sphere/galewsky_bin_output` to see whether your environment is correctly set up.

After this, you could have a look to the file `src/programs/swe_sphere.cpp` which is the starting point for a large variety of benchmarks for the SWE on the sphere.

From hereon, it depends on what you want to do. It's best to get in touch with us (see community information above).


## Unit tests

There's a large set of unit tests available if you execute
``
./run_tests.sh
```
in the root folder.

This tests plenty of features of SWEET, but not everything. It might be necessary to install additional libraries from the local_software folder (see INSTALL file).


## Coding

Please read the coding conventions in doc/sweet/coding_conventions.



## Compilation on cluster:

TODO: Describe how to use MULE
