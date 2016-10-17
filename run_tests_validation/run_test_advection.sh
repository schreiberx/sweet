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
scons --threading=omp --unit-test=test_advection --gui=disable --plane-spectral-space=disable --libfft=disable --mode=release --plane-spectral-dealiasing=disable
EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 1 -N 32 --test-mode 0 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 2 -N 32 --test-mode 0 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 3 -N 32 --test-mode 0 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 4 -N 32 --test-mode 0 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit



echo
echo "***********************************************"
echo "TEST ADVECTION: convergence in time (release) $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_advection --gui=disable --plane-spectral-space=disable --libfft=disable --mode=release --plane-spectral-dealiasing=disable

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 1 -N 32 --test-mode 1 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 2 -N 32 --test-mode 1 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 3 -N 32 --test-mode 1 -G 0 -t 10 -S 0"
echo "$EXEC"
$EXEC || exit

EXEC="./build/test_advection_omp_gnu_release -X 1000000 -Y 1000000 --velocity-u 20000 --velocity-v 0 --advection-scheme 2 --staggered-use-analytical-solution 1 -C 0.1 -H 0 -R 4 -N 32 --test-mode 1 -G 0 -t 10 -S 0"
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

