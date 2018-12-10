#
# MULE Compiler environment
#
export CC=gcc-8
export CXX=g++-8
export F90=gfortran-8
export LINK=$CXX
#export LD=ld

export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpif90

export MULE_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=stdc++

