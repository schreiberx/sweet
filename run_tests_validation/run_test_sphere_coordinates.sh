#! /bin/bash


echo "***********************************************"
echo "Running tests for sphere coordinates"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --unit-test=test_sphere_coordinates"
#SCONS="scons --threading=omp --unit-test=test_sphere_advection --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=release"
echo "$SCONS"
$SCONS || exit 1


EXEC="./build/test_sphere_coordinates*"
echo "$EXEC"
$EXEC || exit
echo
echo
echo



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
