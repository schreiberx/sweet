#! /bin/bash


echo "***********************************************"
echo "Running tests for Operators2d operations"
echo "***********************************************"

cd ..
make clean
scons  --unit-test=test_operators2d --spectral-space=enable --mode=release --spectral-dealiasing=disable
EXEC="./build/test_operators2d_spectral_libfft_gnu_release -N 16 -S 1"
echo "$EXEC"
$EXEC || exit
EXEC="./build/test_operators2d_spectral_libfft_gnu_release -N 16 -S 0"
echo "$EXEC"
$EXEC || exit


make clean
scons  --unit-test=test_operators2d --spectral-space=enable --mode=release --spectral-dealiasing=enable

# test spectral derivatives
EXEC="./build/test_operators2d_spectral_libfft_dealiasing_gnu_release -N 16 -S 1"
echo "$EXEC"
$EXEC || exit

# test cartesian derivatives
EXEC="./build/test_operators2d_spectral_libfft_dealiasing_gnu_release -N 16 -S 0"
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
