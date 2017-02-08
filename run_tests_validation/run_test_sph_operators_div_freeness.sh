#! /bin/bash


echo "**********************************************************"
echo "Running tests for SPH DIV freeness"
echo "This test exists because of issues with the complex-complex div freeness"
echo "**********************************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_sph_operators_div_freeness --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS

./build/test_sph_operators_div_freeness*_debug -M 256 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
