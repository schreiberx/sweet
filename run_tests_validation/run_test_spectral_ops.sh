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

cd ../

X=$MAX_SCALE
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
SCONS="scons --threading=omp --unit-test=test_spectral_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=enable"
$SCONS
./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0 || exit
./build/test_spectral_ops_*_release -n 32 -m 32 -X $X -Y $X -S 1 || exit

X=$MIN_SCALE
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
SCONS="scons --threading=omp --unit-test=test_spectral_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=enable"
$SCONS
./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0 || exit
./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 1 || exit



X=$MAX_SCALE
Y=$MIN_SCALE
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) $X x $Y"
echo "***********************************************"
make clean
SCONS="scons --threading=omp --unit-test=test_spectral_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=disable"
$SCONS
EXEC="./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $Y -S 0"
echo "$EXEC"
$EXEC || exit
EXEC="./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $Y -S 1"
echo "$EXEC"
$EXEC || exit

X=$MAX_SCALE
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) $X"
echo "***********************************************"
make clean
SCONS="scons --threading=omp --unit-test=test_spectral_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=disable"
$SCONS
EXEC="./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 1"
echo "$EXEC"
$EXEC || exit
EXEC="./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0"
$EXEC || exit

X=$MIN_SCALE
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) $X"
echo "***********************************************"
make clean
SCONS="scons --threading=omp --unit-test=test_spectral_ops --gui=disable --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=disable"
$SCONS
./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 1 || exit
./build/test_spectral_ops_*_release -n 128 -m 128 -X $X -Y $X -S 0 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
