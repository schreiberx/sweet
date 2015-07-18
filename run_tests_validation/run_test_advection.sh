#! /bin/bash


echo "***********************************************"
echo "Running tests for Spectral operations"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


# 10^4 km
MAX_SCALE=$((10000*1000))
MIN_SCALE=0.01

cd ..

echo "MAX/MIN SCALE: $MAX_SCALE / $MIN_SCALE"

X=$MAX_SCALE
Y=$MAX_SCALE
echo
echo "***********************************************"
echo "TEST ADVECTION (release) $X"
echo "***********************************************"
make clean
scons --compile-program=test_advection --gui=disable --spectral-space=disable --mode=release --spectral-dealiasing=disable
EXEC="./build/example_test_advection_gnu_release -X $X -Y $Y -N 32 -a 10 -v 0  -C 0.005 -S 0 -R 1 -c 1 -H 0"
echo "$EXEC"
$EXEC || exit
EXEC="./build/example_test_advection_gnu_release -X $X -Y $Y -N 32 -a 10 -v 0  -C 0.005 -S 0 -R 2 -c 1 -H 0"
echo "$EXEC"
$EXEC || exit
EXEC="./build/example_test_advection_gnu_release -X $X -Y $Y -N 32 -a 10 -v 0  -C 0.005 -S 0 -R 3 -c 1 -H 0"
echo "$EXEC"
$EXEC || exit
EXEC="./build/example_test_advection_gnu_release -X $X -Y $Y -N 32 -a 10 -v 0  -C 0.005 -S 0 -R 4 -c 1 -H 0"
echo "$EXEC"
$EXEC || exit




echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
