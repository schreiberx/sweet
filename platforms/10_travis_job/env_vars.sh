
#
# Default (fallback) configuration file
#

echo_warning_hline
echo_warning "This is the environment file for Travis-CI service"
echo_warning "We use Travis to test SWEET"
echo_warning_hline

#
# Compiler environment
#
#export SWEET_F90=gfortran-8
#export SWEET_CC=gcc-8
#export SWEET_CXX=g++-8

# Other variables are automatically set via

# OpenMPI
export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpic++
export SWEET_MPIF90=mpif90
export SWEET_MPILINK=mpic++

#export SWEET_LINK=$SWEET_CXX
#export SWEET_MPILINK=$SWEET_MPICXX

#
# local software compile overrides
#
#export F90=$SWEET_F90
#export FC=$F90
#export CC=$SWEET_CC
#export CXX=$SWEET_CXX
#export LINK=$SWEET_CXX
#export LD=ld

