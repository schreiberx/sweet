#
# How to express dependency in job scheduler
#

export MULE_JOB_SCHEDULER_DEPENDENCY="--dependency=afterany:%JOBID%"

# Limited number of jobs in queue
export MULE_JOB_SCHEDULER_NUM_JOB_LIMITATION=50

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

export CC=gcc
export CXX=g++
export F90=gfortran
export FC=gfortran
export LINK=$CXX
export LD=ld

export MULE_MPICC=mpigcc
export MULE_MPICXX=mpigxx
export MULE_MPIF90=mpifc
export MULE_MPILINK=$MULE_MPICXX

