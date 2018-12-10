
#
# Default (fallback) configuration file
#

echo_warning_hline
echo_warning "This is the environment file for Travis-CI service"
echo_warning "We use Travis to test MULE"
echo_warning_hline

#
# Compiler environment
#
#export F90=gfortran-8
#export CC=gcc-8
#export CXX=g++-8

# Other variables are automatically set via

# OpenMPI
export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpif90

export MULE_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=stdc++
