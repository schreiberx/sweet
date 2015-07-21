#! /bin/bash


echo "***********************************************"
echo "Running tests for Spectral operations"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


cd ..

echo
echo "***********************************************"
echo "TEST ADVECTION: convergence in space (release) $X"
echo "***********************************************"
make clean
scons --unit-test=test_advection --gui=disable --spectral-space=disable --mode=release --spectral-dealiasing=disable
EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 1 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 2 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 3 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 4 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit



echo
echo "***********************************************"
echo "TEST ADVECTION: convergence in time (release) $X"
echo "***********************************************"
make clean
scons --unit-test=test_advection --gui=disable --spectral-space=disable --mode=release --spectral-dealiasing=disable

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 1 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 2 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 3 -N 64 -e 0 -G 0 -t 10"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_gnu_release -a 20000 -c 2 -d 1 -C 0.1 -H 0 -R 4 -N 64 -e 0 -G 0 -t 10"
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
