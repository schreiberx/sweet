
#
# Default GNU configuration file
#

#
# Compiler environment
#
export F90=ifort
export CC=icc
export CXX=icpc
export FC=$F90
export LD=ld
export MULE_LINK=$MULE_CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpicxx
export MULE_MPIF90=mpif90

export MULE_MPILINK=mpicc
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=gfortran

# Load more recent compiler
module unload gcc

#module load gcc/9.3.0
echo "Loading gcc/8.4.0"
module load gcc/8.4.0

# Load this Python module to provide Python.h header file
module load python/3.6_intel

