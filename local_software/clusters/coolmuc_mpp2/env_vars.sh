
#
# Configuration file for CoolMUC mpp2 login nodes
#

echo "Loading SWEET environment specifically for CoolMUC mpp2-* login nodes"


echo "Loading GCC/8"
module unload gcc
module load gcc/8

#echo "Loading binutils"
#module load binutils/2.25

#module unload intel
#module load intel/18.0


# Helper environment variables
# MPIF90 is used by e.g. libpfasst
export SWEET_F90=gfortran
export SWEET_CC=gcc
export SWEET_CXX=g++

export SWEET_MPIF90=mpifc
export SWEET_MPICC=mpigcc
export SWEET_MPICXX=mpigxx

