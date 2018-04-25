#! /bin/bash


echo "***********************************************"
echo "Running tests for Spectral operations"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


# 10^4 km
MAX_SCALE=$((10000*1000))
MIN_SCALE=0.01


echo "MAX/MIN SCALE: $MAX_SCALE / $MIN_SCALE"

cd ..
if false; then
	X=$MAX_SCALE
	Y=$MIN_SCALE
	echo
	echo "***********************************************"
	echo "TEST CART DIFF OPS (release) $X x $Y"
	echo "***********************************************"
	make clean
	scons --threading=omp --unit-test=test_plane_diff_stencil_ops --gui=disable --plane-spectral-space=disable --mode=release --plane-spectral-dealiasing=disable || exit 1
	EXEC="./build/test_plane_diff_stencil_ops_*_release  -X $X -Y $Y -S 0"
	echo "$EXEC"
	$EXEC || exit 1
fi

X=$MAX_SCALE
Y=$MAX_SCALE
echo
echo "***********************************************"
echo "TEST CART DIFF OPS (release) $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_plane_diff_stencil_ops --gui=disable --plane-spectral-space=disable --mode=release --plane-spectral-dealiasing=disable || exit 1
EXEC="./build/test_plane_diff_stencil_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0"
$EXEC || exit 1

X=$MIN_SCALE
echo
echo "***********************************************"
echo "TEST CART DIFF OPS (release) $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_plane_diff_stencil_ops --gui=disable --plane-spectral-space=disable --mode=release --plane-spectral-dealiasing=disable || exit 1
./build/test_plane_diff_stencil_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0 || exit 1

X=$MAX_SCALE
echo
echo "***********************************************"
echo "TEST CART DIFF OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_plane_diff_stencil_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=enable || exit 1
./build/test_plane_diff_stencil_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0 || exit 1

X=$MIN_SCALE
echo
echo "***********************************************"
echo "TEST CART DIFF OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_plane_diff_stencil_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=enable || exit 1
./build/test_plane_diff_stencil_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0 || exit 1



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
