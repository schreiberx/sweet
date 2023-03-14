
#
# Default (fallback) configuration file
#

echo_warning_hline
echo_warning "This is the environment file for GitLab-CI service"
echo_warning "We use GitLab to test MULE and this software development"
echo_warning_hline

# OpenMPI
export MULE_MPICC=mpicc
export MULE_MPICXX=mpic++
export MULE_MPIF90=mpif90

if true; then
	# If we link with mpif90, we have to add stdc++ for C++
	export MULE_MPILINK=mpif90
	export MULE_MPILIBS=stdc++
else
	export MULE_MPILINK=mpic++
	export MULE_MPILIBS=
fi

