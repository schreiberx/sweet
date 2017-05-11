#! /bin/bash


cd ..

echo
echo "***********************************************"
echo "TEST QUADRATURE"
echo "***********************************************"
make clean
scons --unit-test=test_quadrature
EXEC="./build/test_quadrature_*_release"
echo "$EXEC"
$EXEC || exit

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

