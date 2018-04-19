#! /bin/bash


echo "***********************************************"
echo "Running tests for SPH solver"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_sph_solver --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit 1


./build/test_sph_solver*_debug -M 128 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
