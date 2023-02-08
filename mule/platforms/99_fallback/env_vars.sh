
#
# Default (fallback) configuration file
#

echo_warning_hline
echo_warning "This is the fallback environment using GNU compiler environment"
echo_warning "In the case that you just want to test MULE, that's fine."
echo_warning "If you want to actively develop in MULE, think about setting"
echo_warning "up your own platform."
echo_warning_hline

#
# Compiler environment
#
#export F90=gfortran-8
#export CC=gcc-8
#export CXX=g++-8
#export LINK=$CXX
#export LD=ld

export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpif90

export MULE_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=stdc++

export MULE_CC_COMPILER=gcc
export MULE_CXX_COMPILER=gcc
export MULE_F90_COMPILER=gcc
