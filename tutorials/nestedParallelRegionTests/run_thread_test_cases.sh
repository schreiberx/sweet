#! /bin/bash

hline() {
	echo "********************************************************************************"
}

#EXEC=../../build/thread_affinities_COMP_plspec_pldeal_fft_thomp_release
EXEC=../../build/thread_affinities_COMP_thomp_release

run_test() {
	echo

	export OMP_PLACES="$1"
	export OMP_PROC_BIND="$2"

	hline
	echo "OMP_PLACES: $OMP_PLACES"
	echo "OMP_PROC_BIND: $OMP_PROC_BIND"
	hline
	$EXEC
	hline
}


export OMP_MAX_ACTIVE_LEVELS=2
#export OMP_THREADING=4,2


#
# Places:
#
# * socket: threads of a socket
# * cores: sets of all threads of one core
# P="cores"
# * threads: all individual threads
# P="threads"
#


P="{0},{1},{2},{3},{4},{5},{6},{7}"
# Same as before:

NPROC=$(nproc)

NPROC=$((NPROC/2))
P="{0}:$NPROC"

export OMP_NUM_THREADS=$NPROC

#run_test threads spread
#run_test threads close
run_test $P spread,close
run_test $P close,spread


echo "WARNING: Nested parallelization doesn't work with the GNU version. USE LLVM!"
echo "WARNING: Nested parallelization doesn't work with the GNU version. USE LLVM!"
echo "WARNING: Nested parallelization doesn't work with the GNU version. USE LLVM!"

