

BASEDIR=`pwd`
cd ../../

STDEXEC="./build/swe_sph_and_rexi_spherespectral_spheredealiasing_omp_libfft_gnu_release -M 128 -s 9 --nonlinear=0 -g 1 -f 1 -a 1 -H 1 -o 0.5 -t 7 --rexi-sphere-preallocation 0"


if false; then
	for o in 1 2 4; do
		for t in 0.0100 0.0010 0.0001; do
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
		for t in 0.0100 0.0010 0.0001; do
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
		for t in 1.0000 0.1000 0.0100 0.0010; do
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
		for t in 1.0000 0.1000 0.0100 0.0010 0.0001; do
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
	for m in 16 32 64 128 256 512 1024; do
		for t in 0.5 1.0 2.0 3.5 7.0; do
			M=$(printf "%05d" $m)
			OUTDIR="$BASEDIR/output_rexi$M""_ts$t"

			mkdir -p "$OUTDIR" 
			EXEC="$STDEXEC --timestepping-method=100 --rexi-m=$M -C -$t "
			echo "$EXEC"
			rm -f prog_*
			$EXEC || exit 1
			mv prog* "$OUTDIR" || exit 1
		done
	done
fi
