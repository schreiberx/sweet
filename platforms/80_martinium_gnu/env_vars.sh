#
# SWEET Compiler environment
#
export CC=gcc-8
export CXX=g++-8
export F90=gfortran-8
export LINK=$CXX
#export LD=ld

export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpic++
export SWEET_MPIF90=mpif90

export SWEET_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export SWEET_MPILIBS=stdc++

