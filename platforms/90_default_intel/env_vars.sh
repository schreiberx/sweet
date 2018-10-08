
#
# Default GNU configuration file
#


#
# Compiler environment
#
export SWEET_F90=ifort
export SWEET_CC=icc
export SWEET_CXX=icpc

#export SWEET_MPICC=mpicc
#export SWEET_MPICXX=mpic+_
#export SWEET_MPIF90=mpif90

export SWEET_LINK=$SWEET_CXX
#export SWEET_MPILINK=$SWEET_MPICXX

#
# local software compile overrides
#
export F90=$SWEET_F90
export FC=$F90
export CC=$SWEET_CC
export CXX=$SWEET_CXX
export LINK=$SWEET_CXX
export LD=ld

