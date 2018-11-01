
#
# Default (fallback) configuration file
#

echo_warning_hline
echo_warning "This is the fallback environment using GNU compiler environment"
echo_warning "In the case that you just want to test SWEET, that's fine."
echo_warning "If you want to actively develop in SWEET, think about setting"
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

export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpic++
export SWEET_MPIF90=mpif90

export SWEET_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export SWEET_MPILIBS=stdc++
