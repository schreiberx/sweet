#! /bin/bash


echo "***********************************************"
echo "Running tests for SPH operators"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_sph_operators --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit 1

time ./build/test_sph_operators*_debug -M 256 --reuse-plans=-1 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
