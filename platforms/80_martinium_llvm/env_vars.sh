#
# SWEET Compiler environment
#
export SWEET_CC=clang-7
export SWEET_CXX=clang++-7
export SWEET_F90=

export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpic++
export SWEET_MPIF90=mpif90

export SWEET_LINK=$SWEET_CXX
export SWEET_MPILINK=$SWEET_MPICXX

#
# local software compile overrides
#
export F90=$SWEET_F90
export FC=$F90
export CC=$SWEET_CC
export CXX=$SWEET_CXX
export LINK=$SWEET_CXX
export LD=ld

