#! /bin/bash


echo "***********************************************"
echo "Running tests for interpolation on the plane"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_plane_interpolation --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS

./build/test_plane_interpolation*_debug -M 32 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
