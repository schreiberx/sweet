#! /bin/bash


echo "***********************************************"
echo "Running tests for sampler operations"
echo "***********************************************"

cd ..

echo
echo "***********************************************"
echo "TEST SPECTRAL OPS (release) ALIASING CONTROL $X"
echo "***********************************************"
make clean
scons --threading=omp --unit-test=test_samplers --gui=disable --plane-spectral-dealiasing=disable
./build/test_samplers_libfft_omp_gnu_release -s 2 -n 128 -m 128 || exit
./build/test_samplers_libfft_omp_gnu_release -s 3 -n 128 -m 128 || exit



echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
