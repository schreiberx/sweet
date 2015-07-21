#! /bin/bash


echo "***********************************************"
echo "Running tests for Spectral operations"
echo "***********************************************"

export COMPILER=gnu
export OMP_PROC_BIND=close


# 10^4 km
MAX_SCALE=$((10000*1000))

MIN_SCALE=0.01

echo "MAX/MIN SCALE: $MAX_SCALE / $MIN_SCALE"

cd ..

# earth-like
X=$((40000*1000))
Y=$((20000*1000))
SCALE=1000
RES_X=$((40*$SCALE))
RES_Y=$((20*$SCALE))
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
scons --compiler=$COMPILER --unit-test=test_spectral_ops --gui=disable --spectral-space=enable --mode=release --spectral-dealiasing=enable
./build/test_spectral_ops_spectral_dealiasing_gnu_release -n $RES_X -m $RES_Y -X $X -Y $Y -S 0 || exit
./build/test_spectral_ops_spectral_dealiasing_gnu_release -n $RES_X -m $RES_Y -X $X -Y $Y -S 1 || exit


# unit-sphere-like
X=2
Y=1
SCALE=1000
RES_X=$((40*$SCALE))
RES_Y=$((20*$SCALE))
echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
scons --compiler=$COMPILER --unit-test=test_spectral_ops --gui=disable --spectral-space=enable --mode=release --spectral-dealiasing=enable
./build/test_spectral_ops_spectral_dealiasing_gnu_release -n $RES_X -m $RES_Y -X $X -Y $Y -S 0 || exit
./build/test_spectral_ops_spectral_dealiasing_gnu_release -n $RES_X -m $RES_Y -X $X -Y $Y -S 1 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
