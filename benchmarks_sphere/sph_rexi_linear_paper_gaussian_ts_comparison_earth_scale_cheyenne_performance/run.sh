#! /bin/bash

BASEDIR=`pwd`
SWEETDIR="`pwd`/../../"


# simtime: 7
# output: 14 times
# ts: about 0.01

# simtime: 1.5 days (129600 seconds)
# output: 12 times (10800 seconds)
# ts: about 200

OUTPUT_TIME=10800
SIMTIME=129600

STDEXEC="$SWEETDIR/build/swe_sphere_spherespectral_spheredealiasing_omp_libfft_gnu_release -M 128 -s 9 --nonlinear=0 -o $OUTPUT_TIME -t $SIMTIME --rexi-sphere-preallocation 0 --use-robert-functions=1 -F 0"

# RK for video
if true; then
#if false; then
	for o in 4; do
		for t in 050; do
			cd "$BASEDIR"
			OUTDIR="$BASEDIR/output_video_rk$o""_ts$t"
			mkdir -p "$OUTDIR" 
			cd "$OUTDIR"

			EXEC="$STDEXEC --timestepping-method=l_erk --timestepping-order=$o -C -$t -o 100"
			echo "$EXEC"
			rm -f output_prog_*
			$EXEC || exit 1
		done
	done
fi

exit 1

# RK
#if true; then
if false; then
	for o in 1 2 4; do
		for t in 050 100 200 400 800; do
			cd "$BASEDIR"
			OUTDIR="$BASEDIR/output_rk$o""_ts$t"
			mkdir -p "$OUTDIR" 
			cd "$OUTDIR"

			EXEC="$STDEXEC --timestepping-method=l_erk --timestepping-order=$o -C -$t" # -o $((10800/100))"
			echo "$EXEC"
			rm -f output_prog_*
			$EXEC || exit 1
		done
	done
fi

# Leapfrog
#if true; then
if false; then
	for o in 2; do
		for t in 050 100 200 400 800; do
			cd "$BASEDIR"
			OUTDIR="$BASEDIR/output_lf$o""_ts$t"
			mkdir -p "$OUTDIR" 
			cd "$OUTDIR"

			EXEC="$STDEXEC --timestepping-method=l_lf --timestepping-order=$o -C -$t"
			echo "$EXEC"
			rm -f output_prog_*
			$EXEC || exit 1
		done
	done
fi


# Euler implicit
#if true; then
if false; then
	for o in 1; do
		for t in 050 100 200 400 800; do
			cd "$BASEDIR"
			OUTDIR="$BASEDIR/output_cn2""_ts$t"
			mkdir -p "$OUTDIR" 
			cd "$OUTDIR"

			EXEC="$STDEXEC --timestepping-method=l_cn --timestepping-order=2 -C -$t"
			echo "$EXEC"
			rm -f output_prog_*
			$EXEC || exit 1
		done
	done
fi



# REXI
if true; then
	for m in 16 32 512; do
		for t in 800; do
			cd "$BASEDIR"
			M=$(printf "%05d" $m)
			T=$(printf "%06d" $t)
			OUTDIR="$BASEDIR/output_rexi$M""_rexinorm$n""_ts$T"
			mkdir -p "$OUTDIR" 
			cd "$OUTDIR"

			EXEC="$STDEXEC --timestepping-method=l_rexi --rexi-m=$M --rexi-normalization=$n -C -$t "
			echo "$EXEC"
			rm -f output_prog_*
			$EXEC || exit 1
		done
	done
fi


# REXI
if true; then
	for m in 1024 4192; do
		for t in 129600; do
			cd "$BASEDIR"
			M=$(printf "%05d" $m)
			T=$(printf "%06d" $t)
			OUTDIR="$BASEDIR/output_rexi$M""_rexinorm$n""_ts$T"
			mkdir -p "$OUTDIR" 
			cd "$OUTDIR"

			EXEC="$STDEXEC --timestepping-method=l_rexi --rexi-m=$M --rexi-normalization=$n -C -$t "
			echo "$EXEC"
			rm -f output_prog_*
			$EXEC || exit 1
		done
	done
fi

