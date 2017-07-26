#! /bin/bash


echo "***********************************************"
echo "Running tests for plane restriction and interpolation"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ..

make clean
SCONS="scons --threading=omp --unit-test=test_plane_restriction_interpolation --mode=debug --plane-spectral-dealiasing=disable"
$SCONS

./build/test_plane_restriction_interpolation_planespectral_omp_libfft_gnu_debug  -N 64 || exit 1



make clean
SCONS="scons --threading=omp --unit-test=test_plane_restriction_interpolation --mode=debug --plane-spectral-dealiasing=enable"
$SCONS

./build/test_plane_restriction_interpolation_planespectral_planedealiasing_omp_libfft_gnu_debug  -N 64 || exit 1


echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
