#! /bin/bash


echo "***********************************************"
echo "Running tests for SPH solver complex"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_sph_solver_complex --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=debug --quadmath=enable"
echo "$SCONS"
$SCONS || exit 1


for r in 0.1 10000; do

	EXEC="./build/test_sph_solver_complex*_debug -f 1 -g 1 -H 1 -a $r -M 256"
	echo "$EXEC"
	$EXEC || exit

done



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
