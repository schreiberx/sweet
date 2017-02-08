

BASEDIR=`pwd`
cd ../../

STDEXEC="./build/swe_sph_and_rexi_spherespectral_spheredealiasing_omp_libfft_gnu_release -M 128 --timestepping-order=4 -s 6 --nonlinear=0 -g 1 -f 1 -a 1 -H 1 -o 0.5 -t 7 --rexi-sphere-preallocation 0"


if false; then
	OUTDIR="$BASEDIR/output_rk4"
	mkdir -p "$OUTDIR" 
	EXEC="$STDEXEC --timestepping-method=1 -C -0.01"
	echo "$EXEC"
	rm -f prog_*
	$EXEC || exit 1
	mv prog* "$OUTDIR" || exit 1
fi


#for m in 16 32 64 128 256 512 1024; do
for m in 1024 2048; do

	M=$(printf "%05d" $m)
	OUTDIR="$BASEDIR/output_rexi$M"

	mkdir -p "$OUTDIR" 
	EXEC="$STDEXEC --timestepping-method=100 --rexi-m=$M -C -0.5 "
	echo "$EXEC"
	rm -f prog_*
	$EXEC || exit 1
	mv prog* "$OUTDIR" || exit 1
done
