#! /bin/bash


echo "***********************************************"
echo "Running tests for Spectral operations"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST RK TIMESTEPPING (release) $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_runge_kutta --gui=disable --spectral-space=disable --libfft=disable --mode=release --spectral-dealiasing=disable
EXEC="./build/test_runge_kutta_omp_gnu_release -N 128 -G 0 -C -1 -t 5 --function-order 4 -S 0"
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
