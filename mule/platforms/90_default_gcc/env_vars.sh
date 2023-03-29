
#
# Default GNU configuration file
#

#
# Compiler environment
#

if [ "`uname`" == "Darwin" ]; then
        # Setup environment for default compiler (installed via `brew install gcc@12`)
        export PYTHON=python3
        export CC=gcc-12
        export CXX=g++-12
        export F90=gfortran-12
else
	test -z "$F90" && export F90=gfortran
	test -z "$CC" && export CC=gcc
	test -z "$CXX" && export CXX=g++
	test -z "$FC" && export FC=$F90
	test -z "$LD" && export LD=ld
fi


export MULE_LINK=$MULE_CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpif90

if true; then
	# If we link with mpif90, we have to add stdc++ for C++
	export MULE_MPILINK=mpif90
	export MULE_MPILIBS=stdc++
else
	export MULE_MPILINK=mpic++
	export MULE_MPILIBS=
fi


export MULE_CC_COMPILER=gcc
export MULE_CXX_COMPILER=gcc
export MULE_F90_COMPILER=gcc

