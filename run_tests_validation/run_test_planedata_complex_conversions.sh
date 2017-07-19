#! /bin/bash


echo "***********************************************"
echo "Running tests for plane data conversions"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_planedata_complex_conversions --gui=disable --plane-spectral-space=enable --plane-spectral-dealiasing=disable --mode=debug"
echo "$SCONS"
$SCONS

./build/test_planedata_complex_conversions_*debug -M 16 || exit




make clean
SCONS="scons --threading=omp --unit-test=test_planedata_complex_conversions --gui=disable --plane-spectral-space=enable --plane-spectral-dealiasing=enable --mode=debug"
echo "$SCONS"
$SCONS

./build/test_planedata_complex_conversions_*debug -M 16 || exit


echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
