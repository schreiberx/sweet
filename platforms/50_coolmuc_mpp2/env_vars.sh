
#
# Configuration file for CoolMUC mpp2 login nodes
#

echo "Loading GCC/8"
module unload gcc
module load gcc/8

#echo "Loading binutils"
#module load binutils/2.25

#module unload intel
#module load intel/18.0


#
# Compiler environment
#
export SWEET_F90=gfortran
export SWEET_CC=gcc
export SWEET_CXX=g++

export SWEET_MPICC=mpigcc
export SWEET_MPICXX=mpigxx
export SWEET_MPIF90=mpifc

