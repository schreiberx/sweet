#
# Default configuration file
#



# Load more recent compiler

module load gcc

#
# Don't touch intel stuff, otherwise the mpicxx is not available
#
module load intel


#
# Load this Python module to provide Python.h header file
# The regular python3 include directory is messed up
#
module unload python
module load python/3.6_intel


export MULE_JOB_SCHEDULER_NUM_JOB_LIMITATION=30
export MULE_JOB_SCHEDULER_SLEEP_SECS=30


#
# Compiler environment
#
export F90=ifort
export CC=icc
export CXX=icpc
export FC=$F90
export LD=ld
export MULE_LINK=$MULE_CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpicxx
export MULE_MPIF90=mpif90


export MULE_MPILINK=mpif90
# If we link with mpif90, we have to add stdc++ for C++
export MULE_MPILIBS=stdc++
#export MULE_MPILIBS=gfortran


export MULE_CC_COMPILER=intel
export MULE_CXX_COMPILER=intel
export MULE_F90_COMPILER=intel

export MULE_USER_EMAIL=
