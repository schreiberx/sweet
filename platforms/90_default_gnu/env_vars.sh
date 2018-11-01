
#
# Default GNU configuration file
#

#
# Compiler environment
#
export F90=gfortran
export CC=gcc
export CXX=g++
export FC=$F90
export LD=ld
export SWEET_LINK=$SWEET_CXX

export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpic++
export SWEET_MPIF90=mpif90

export SWEET_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export SWEET_MPILIBS=stdc++
