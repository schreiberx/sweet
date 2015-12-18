#! /bin/bash

CURDIR=`pwd`

cd ../../

. ./local_software/env_vars.sh
scons --compiler=gnu --program=swe_rexi --spectral-space=disable --libfft=enable --rexi-parallel-sum=disable --spectral-dealiasing=disable --mode=release --threading=omp --numa-block-allocator=1 || exit 1
BIN=../../build/swe_rexi_libfft_omp_numaallocator1_gnu_release

cd "$CURDIR"

SIMTIME=50
OUTPUT_INTERVAL=1
CFL=0.3

RES_LIST="32 64 128 256"


for N in $RES_LIST; do
	echo "WARNING!!!!!! The analytical results are computed on an A-grid - only compare the height!"
	echo "Generating analytical correct solution with resolution $N x $N"
	$BIN -N $N -X 1 -Y 1 -g 1 -H 1 -f 1 -t $SIMTIME --compute-error=1 --timestepping-mode=2 --use-fd-for-complex-array=1 -S 0 -s 5 -v 2 -O output_solution_$N -o $OUTPUT_INTERVAL -C -$OUTPUT_INTERVAL > "output_solution_$N""_diagnostics"

	echo "Running FD simulation with resolution $N x $N"
	$BIN -N $N -X 1 -Y 1 -g 1 -H 1 -f 1 -R 4 -t $SIMTIME --compute-error=1 --staggering=1 --timestepping-mode=0 --use-fd-for-complex-array=1 -S 0 -s 5 -v 2 -O output_$N -o $OUTPUT_INTERVAL -C $CFL > "output_$N""_diagnostics"
done
