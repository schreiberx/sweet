#
# Tag in header of job subscription files to express dependency to another job
# This is highly important for the plan generation of spectral transformations
#
export JOB_SCHEDULER_DEPENDENCY="-W depend=afterany:%JOBID%"


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
#export SWEET_CC=gcc
#export SWEET_CXX=g++
#export SWEET_F90=gfortran

export SWEET_CC=/glade/u/apps/ch/opt/gnu/8.1.0/bin/gcc
export SWEET_CXX=/glade/u/apps/ch/opt/gnu/8.1.0/bin/g++
export SWEET_F90=/glade/u/apps/ch/opt/gnu/8.1.0/bin/gfortran

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


