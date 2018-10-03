#
# Configuration file for CoolMUC mpp2 login nodes
#


#
# Tags in header of batch files
#
# This is important for the SHTNS plan generation scripts
#
export BATCH_FILE_TAG="#SBATCH"


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
export SWEET_CC=gcc
export SWEET_CXX=g++
export SWEET_F90=gfortran

export SWEET_MPICC=mpigcc
export SWEET_MPICXX=mpigxx
export SWEET_MPIF90=mpifc

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

