
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
#export F90=gfortran-8
#export CC=gcc-8
#export CXX=g++-8

# Other variables are automatically set via

# OpenMPI
export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpic++
export SWEET_MPIF90=mpif90
export SWEET_MPILINK=mpic++

