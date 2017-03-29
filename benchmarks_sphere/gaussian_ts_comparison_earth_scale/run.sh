

BASEDIR=`pwd`
cd ../../


# simtime: 7
# output: 14 times
# ts: about 0.01

# simtime: 1.5 days (129600 seconds)
# output: 12 times (10800 seconds)
# ts: about 200

OUTPUT_TIME=10800
SIMTIME=129600

STDEXEC="./build/swe_sph_and_rexi_spherespectral_spheredealiasing_omp_libfft_gnu_release -M 128 -s 9 --nonlinear=0 -o $OUTPUT_TIME -t $SIMTIME --rexi-sphere-preallocation 0 --use-robert-functions=1 -F 0"

# RK
#if true; then
if false; then
	for o in 1 2 4; do
		for t in 050 100 200 400 800; do
			for PDEID in 1; do
				OUTDIR="$BASEDIR/output_rk$o""_ts$t""_pde$PDEID"
				mkdir -p "$OUTDIR" 
				EXEC="$STDEXEC --timestepping-method=1 --timestepping-order=$o --pde-id=$PDEID -C -$t" # -o $((10800/100))"
				echo "$EXEC"
				rm -f output_prog_*
				$EXEC || exit 1
				mv output_prog* "$OUTDIR"
			done
		done
	done
fi

# Leapfrog
#if true; then
if false; then
	for o in 2; do
		for t in 050 100 200 400 800; do
			for PDEID in 1; do
				OUTDIR="$BASEDIR/output_lf$o""_ts$t""_pde$PDEID"
				mkdir -p "$OUTDIR" 
				EXEC="$STDEXEC --timestepping-method=2 --timestepping-order=$o --pde-id=$PDEID -C -$t"
				echo "$EXEC"
				rm -f output_prog_*
				$EXEC || exit 1
				mv output_prog* "$OUTDIR"
			done
		done
	done
fi


# Euler implicit
#if true; then
if false; then
	#for o in 1 2; do
	for o in 2; do
		for PDEID in 1; do
			for t in 050 100 200 400 800; do
				OUTDIR="$BASEDIR/output_cn$o""_ts$t""_pde$PDEID"
				mkdir -p "$OUTDIR" 
				EXEC="$STDEXEC --timestepping-method=3 --timestepping-order=$o --pde-id=$PDEID -C -$t"
				echo "$EXEC"
				rm -f output_prog_*
				$EXEC || exit 1
				mv output_prog* "$OUTDIR"
			done
		done
	done
fi



# REXI
#if true; then
if false; then
#	for m in 16 32 64 128 256 512 1024 2048; do
	for m in 16 32 64 128 256 512 1024 2048 4096; do
#	for m in 2048 4096; do
#		for t in 129600 64800 32400 16200 8100 4050 2000 1400 800 400 200 100 050; do
		for t in 129600 64800 32400 16200 8100 4050 2000 1400 800; do
			for PDEID in 1; do
				M=$(printf "%05d" $m)
				T=$(printf "%06d" $t)
				OUTDIR="$BASEDIR/output_rexi$M""_rexinorm$n""_ts$T""_pde$PDEID"

				mkdir -p "$OUTDIR" 
				EXEC="$STDEXEC --timestepping-method=100 --rexi-m=$M --rexi-normalization=$n --pde-id=$PDEID -C -$t "
				echo "$EXEC"
				rm -f output_prog_*
				$EXEC || exit 1
				mv output_prog* "$OUTDIR"
			done
		done
	done
fi


if true; then
#if false; then
	for m in 1024 2048 4192; do
		for t in 800; do
			for PDEID in 1; do
				M=$(printf "%05d" $m)
				T=$(printf "%06d" $t)
				OUTDIR="$BASEDIR/output_rexi$M""_rexinorm$n""_ts$T""_pde$PDEID"

				mkdir -p "$OUTDIR" 
				EXEC="$STDEXEC --timestepping-method=100 --rexi-m=$M --rexi-normalization=$n --pde-id=$PDEID -C -$t "
				echo "$EXEC"
				rm -f output_prog_*
				$EXEC || exit 1
				mv output_prog* "$OUTDIR"
			done
		done
	done
fi

