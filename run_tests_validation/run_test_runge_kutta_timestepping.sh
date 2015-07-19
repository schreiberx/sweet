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
scons --compile-program=test_runge_kutta --gui=enable --spectral-space=disable --mode=release --spectral-dealiasing=disable
EXEC="./build/example_test_runge_kutta_gui_gnu_release -N 32 -G 0 -C -1 -t 5 -a 4"
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
