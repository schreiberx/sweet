#
# Configuration file for CoolMUC mpp2 login nodes
#


#
# Tags in header of batch files
#
# This is important for the SHTNS plan generation scripts
#
export BATCH_FILE_TAG="#PBS"


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


