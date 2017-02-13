

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

STDEXEC="./build/swe_sph_and_rexi_spherespectral_spheredealiasing_omp_libfft_gnu_release -M 128 -s 6 --nonlinear=0 -o $OUTPUT_TIME -t $SIMTIME --rexi-sphere-preallocation 0"


if false; then
	for o in 1 2 4; do
		for t in 050 100 200 400 800; do
			OUTDIR="$BASEDIR/output_rk$o""_ts$t"
			mkdir -p "$OUTDIR" 
			EXEC="$STDEXEC --timestepping-method=1 --timestepping-order=$o -C -$t"
			echo "$EXEC"
			rm -f prog_*
			$EXEC || exit 1
			mv prog* "$OUTDIR" || exit 1
		done
	done
fi


if false; then
	for o in 2; do
		for t in 050 100 200 400; do
			OUTDIR="$BASEDIR/output_lf$o""_ts$t"
			mkdir -p "$OUTDIR" 
			EXEC="$STDEXEC --timestepping-method=2 --timestepping-order=$o -C -$t"
			echo "$EXEC"
			rm -f prog_*
			$EXEC || exit 1
			mv prog* "$OUTDIR" || exit 1
		done
	done
fi


if false; then
	for o in 2; do
		for t in 050 100 200 400; do
			OUTDIR="$BASEDIR/output_cn$o""_ts$t"
			mkdir -p "$OUTDIR" 
			EXEC="$STDEXEC --timestepping-method=3 --timestepping-order=$o -C -$t"
			echo "$EXEC"
			rm -f prog_*
			$EXEC || exit 1
			mv prog* "$OUTDIR" || exit 1
		done
	done
fi


if false; then
	for o in 1; do
		for t in 050 100 200 400; do
			OUTDIR="$BASEDIR/output_irk$o""_ts$t"
			mkdir -p "$OUTDIR" 
			EXEC="$STDEXEC --timestepping-method=3 --timestepping-order=$o -C -$t"
			echo "$EXEC"
			rm -f prog_*
			$EXEC || exit 1
			mv prog* "$OUTDIR" || exit 1
		done
	done
fi


if true; then
	for t in 129600 64800 32400 16200 8100 4050; do
		for m in 16 32 64 128 256 512 1024 2048; do
			for n in 1 0; do
				M=$(printf "%05d" $m)
				T=$(printf "%06d" $t)
				OUTDIR="$BASEDIR/output_rexi$M""_rexinorm$n""_ts$T"

				mkdir -p "$OUTDIR" 
				EXEC="$STDEXEC --timestepping-method=100 --rexi-m=$M --rexi-normalization=$n -C -$t "
				echo "$EXEC"
				rm -f prog_*
				$EXEC || exit 1
				mv prog* "$OUTDIR" || exit 1
			done
		done
	done
fi



