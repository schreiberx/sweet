#
# Tag in header of job subscription files to express dependency to another job
# This is highly important for the plan generation of spectral transformations
#
export MULE_JOB_SCHEDULER_DEPENDENCY="-W depend=afterany:%JOBID%"


MODULES="gnu/8.1.0"
for m in $MODULES; do
	echo
	echo "Loading $m"
	module load $m
done


#
# Compiler environment
#


#
# DO NOT USE icc, icpc or ifort directly on Cheyenne!
# These are python scripts which mess around with the libraries!!!
#
# Use e.g.
#	g++ --show
# to see what's really executed for the compiler
#

export CC=/glade/u/apps/ch/opt/gnu/8.1.0/bin/gcc
export CXX=/glade/u/apps/ch/opt/gnu/8.1.0/bin/g++
export F90=/glade/u/apps/ch/opt/gnu/8.1.0/bin/gfortran
export FC=/glade/u/apps/ch/opt/gnu/8.1.0/bin/gfortran
export LINK=$CXX
export LD=ld

export MULE_MPICC=mpicc
export MULE_MPICXX=mpicxx
export MULE_MPIF90=mpif90
export MULE_MPILINK=$MULE_MPICXX

export MULE_CC_COMPILER=gcc
export MULE_CXX_COMPILER=gcc
export MULE_F90_COMPILER=gcc
