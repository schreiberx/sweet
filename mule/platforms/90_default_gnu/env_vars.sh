
#
# Default GNU configuration file
#

#
# Compiler environment
#

if [ "`uname`" == "Darwin" ]; then
        # Setup environment for default compiler (installed via `brew install gcc@11`)
        export PYTHON=python3
        export CC=gcc-11
        export CXX=g++-11
        export F90=gfortran-11
else
	export F90=gfortran
	export CC=gcc
	export CXX=g++
	export FC=$F90
	export LD=ld
fi

export MULE_LINK=$MULE_CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpif90

export MULE_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=stdc++

