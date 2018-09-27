
#
# Configuration file for CoolMUC mpp2 login nodes
#

if [ "${HOSTNAME:0:9}" == "martinium" ]; then

	echo "Loading SWEET environment specifically for test 'martinium' system"

	#
	# Compiler environment
	#
	export SWEET_F90=gfortran-8
	export SWEET_CC=gcc-8
	export SWEET_CPP=g++-8

	export SWEET_MPICC=mpicc
	export SWEET_MPICXX=mpic+_
	export SWEET_MPIF90=mpif90

fi
