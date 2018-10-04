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
$SCONS || exit 1

EXEC="$(eval echo ./build/test_planedata_complex_conversions_*debug -M 16)"
echo "$EXEC"
$EXEC || exit




make clean
SCONS="scons --threading=omp --unit-test=test_planedata_complex_conversions --gui=disable --plane-spectral-space=enable --plane-spectral-dealiasing=enable --mode=debug"
echo "$SCONS"
$SCONS || exit 1

EXEC="$(eval echo ./build/test_planedata_complex_conversions_*debug -M 16)"
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
