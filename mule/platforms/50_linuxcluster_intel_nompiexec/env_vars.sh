#
# Default GNU configuration file
#


# Load more recent compiler
module unload gcc
module load gcc/8.4.0
#
# Load this Python module to provide Python.h header file
# The regular python3 include directory is messed up
#


# Don't touch intel stuff, otherwise the mpicxx is not availble
#module unload intel
#module load intel-mpi/2019.6.166
#module load intel/19.0.5
#module load intel-mkl/2019.5.281


module unload python
module load python/3.6_intel

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


export MULE_MPILINK=mpicxx
# If we link with mpif90, we have to add stdc++ for C++
#export MULE_MPILIBS=stdc++
export MULE_MPILIBS=gfortran

