#
# Default configuration file
#

module load mpi


#
# Compiler environment
#
export F90=gfortran
export CC=gcc
export CXX=g++
export FC=$F90
export LD=ld
export MULE_LINK=$CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpifort


export MULE_MPILINK=mpicc
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=stdc++

