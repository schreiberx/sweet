#! /bin/bash


NUM_PROCS=`cat /proc/cpuinfo  | grep processor | wc -l`
PROC_RANGE=`seq 0 $NUM_PROCS`

if [ -n "$1" ]; then
	PROC_RANGE=$1
fi

for t in $PROC_RANGE; do
	WISDOM_FILE="FFTW_WISDOM_T$t"
	echo "Creating wisdom in file $WISDOM_FILE"

	if [ -n "$2" ]; then
		WISDOM_FILE="$2"
	fi

	PLANS=""

	# Generate all these plans.
	# cif and cf plans are required to support different platforms.
	# On some platforms, only cif seems to work, on others only cf
#	for n in 8 16 32 64 128; do
	for n in 8 16 32 64 128 256 512 1024; do
		PLANS="$PLANS rf""$n""x""$n"
		PLANS="$PLANS rb""$n""x""$n"
		PLANS="$PLANS cif""$n""x""$n"
		PLANS="$PLANS cib""$n""x""$n"
		PLANS="$PLANS cf""$n""x""$n"
		PLANS="$PLANS cb""$n""x""$n"
	done

	T=""
	if [ $t -ge 1 ]; then 
		T=" -T $t "
		declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
		declare -x OMP_NUM_THREADS=$t
	else
		unset KMP_AFFINITY
		unset OMP_NUM_THREADS
	fi
	FFTW_WISDOM_EXEC="fftw-wisdom $T $PLANS"
	echo "$FFTW_WISDOM_EXEC"

	$FFTW_WISDOM_EXEC > $WISDOM_FILE
done
