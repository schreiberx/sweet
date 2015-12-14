#! /bin/bash


NUM_PROCS=`cat /proc/cpuinfo  | grep processor | wc -l`
for t in `seq 0 $NUM_PROCS`; do
	WISDOM_FILE="FFTW_WISDOM_T$t"
	echo "Creating wisdom in file $WISDOM_FILE"

	PLANS=""
	for n in 8 16 32 64 128 256 512; do
		PLANS="$PLANS rf""$n""x""$n"
		PLANS="$PLANS rb""$n""x""$n"
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
