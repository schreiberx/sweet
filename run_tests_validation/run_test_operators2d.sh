#! /bin/bash


echo "***********************************************"
echo "Running tests for Operators2d operations"
echo "***********************************************"


echo "***********************************************"
echo "Running tests for Operators2d operations > SPECTRAL SPACE, SPECTRAL OPERATORS"
echo "***********************************************"
cd ..
make clean
scons  --unit-test=test_operators2d --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=disable
EXEC="./build/test_operators2d_*_release -N 16 -S 1 --timestepping-method=1 --timestepping-order 4"
echo "$EXEC"
$EXEC || exit

echo "***********************************************"
echo "Running tests for Operators2d operations > SPECTRAL SPACE, STENCIL OPERATORS"
echo "***********************************************"
EXEC="./build/test_operators2d_*_release -N 16 -S 0 --timestepping-method=1 --timestepping-order 4"
echo "$EXEC"
$EXEC || exit


make clean
scons  --unit-test=test_operators2d --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=enable

echo "***********************************************"
echo "Running tests for Operators2d operations > SPECTRAL SPACE, SPECTRAL OPERATORS, DEALIASING"
echo "***********************************************"
# test spectral derivatives
EXEC="./build/test_operators2d_*_release -N 16 -S 1 --timestepping-method=1 --timestepping-order 4"
echo "$EXEC"
$EXEC || exit

# test cartesian derivatives
echo "***********************************************"
echo "Running tests for Operators2d operations > SPECTRAL SPACE, STENCIL OPERATORS, DEALIASING"
echo "***********************************************"
EXEC="./build/test_operators2d_*_release -N 16 -S 0 --timestepping-method=1 --timestepping-order 4"
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
