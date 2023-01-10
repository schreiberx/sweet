
READ THIS CAREFULLY!!! DON'T SKIP ANY PARTS!!!


# Installation information about SWEET

There are scripts in the directory
```
local_software
```
to make handling external software requirements easier.

Not all software packages are required on all platforms. E.g. installing
`autoconf` is only required if the version on the platform is too old.


## First steps

Change the directory to SWEET's root directory, e.g.
```
$ cd sweet
```


### 1. Setup environment

Before doing anything, you *must* setup your environment variable via

```
$ source ./activate.sh
```

### 2. Miniconda

DEACTIVATE ALL (ANA)CONDA ENVIRONMENTS!

Then, install ANACONDA by typing
```
$ cd local_software
$ ./install_miniconda.sh
```
This is not always required, but you're on the safer side to have a Python version and its modules in a particular version.


### 3. 3rd party libraries

Within the 3rd party software directory in local_software,
there are installation scripts named as follows:
```
   ./install_*.sh
```
These are scripts which build and install the 3rd party software automatically.

See next section for the recommended packages


## Recommended packages

This is a list of recommended scripts to trigger for installation.
First of all, include the environment variables which make the
installed software available to the other install scripts:

Once going into the local_software directory
```
$ cd ./local_software
```
install the following packages if required:


Use pip to install other packages:
```
$ pip3 install matplotlib sympy mpmath
```

If your compiler is older than gcc 5.3 (check with gcc --version), then
```
$ ./install_gcc8.2.sh	# GNU compiler
```

The other librar(ies)/y can be installed via:
```
$ ./install_scons.sh	# Makefile replacement
$ ./install_numa.sh	# NUMA aware memory allocation
```

Packages required for simulations on the plane:
```
$ ./install_fftw3.sh	# Fast Fourier library
```

After installing the above mentioned software, you should
be able to compile an example program such as with

```
   $ scons --program=swe_plane
```


Packages required for simulations on the sphere:
```
$ ./install_shtns.sh	# Spherical Harmonics library
$ ./install_lapack.sh # Linear algebra package
```

Packages required for pfasst++ (C++ PFASST):
```
$ ./install_eigen3.sh	# Eigen library
```

Packages required for libpfasst (Fortran PFASST):
```
$ ./install_libpfasst.sh
```

If there are any compilation problems, please send the last
lines of the output and more information on the used system to

Martin Schreiber <schreiberx@gmail.com>


## What if something goes wrong?

Get in touch with us (see community information in the README.md file or drop an
email to schreiberx@gmail.com )


## Mac users

Check out the INSTALL_MACOSX.md file
