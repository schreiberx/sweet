#
# Configuration file for CoolMUC mpp2 login nodes
#


#
# Tags in header of batch files
#
# This is important for the SHTNS plan generation scripts
#
export BATCH_FILE_TAG="#SBATCH"


MODULES="gnu/8.1.0"
for m in $MODULES; do
	echo
	echo "Loading $m"
	module load $m
done


#
# Compiler environment
#

export SWEET_CC=gcc
export SWEET_CXX=g++
export SWEET_F90=gfortran

export SWEET_MPICC=mpicc
export SWEET_MPICXX=mpicxx
export SWEET_MPIF90=mpif90

