#
# Tag in header of job subscription files to express dependency to another job
# This is highly important for the plan generation of spectral transformations
#
export MULE_JOB_SCHEDULER_DEPENDENCY="-W depend=afterany:%JOBID%"


MODULES="intel/18.0.1"
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
# icpc --show
# to see what's really executed for the compiler
#
export CC=/glade/u/apps/opt/intel/2018u1/compilers_and_libraries/linux/bin/intel64/icc
export CXX=/glade/u/apps/opt/intel/2018u1/compilers_and_libraries/linux/bin/intel64/icpc
export F90=/glade/u/apps/opt/intel/2018u1/compilers_and_libraries/linux/bin/intel64/ifort
export FC=$F90
export LINK=$CXX
export LD=ld

export MULE_MPICC=/opt/sgi/mpt/mpt-2.15/bin/mpicc
export MULE_MPICXX=/opt/sgi/mpt/mpt-2.15/bin/mpicxx
export MULE_MPIF90=/opt/sgi/mpt/mpt-2.15/bin/mpif90
export MULE_MPILINK=$MULE_MPICXX

export MULE_CC_COMPILER=intel
export MULE_CXX_COMPILER=intel
export MULE_F90_COMPILER=intel
