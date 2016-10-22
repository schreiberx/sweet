#! /bin/bash


echo "***********************************************"
echo "Running tests for antialiasing patterns"
echo "***********************************************"

cd ../


if false; then

	make clean

	SCONS="scons --unit-test=test_antialiasing_patterns --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=disable"
	echo "$SCONS"
	$SCONS
	EXEC="./build/test_antialiasing_patterns_*_release -N 16 -S 1"
	echo "$EXEC"
	$EXEC || exit

	EXEC="./build/test_antialiasing_patterns_*_release -N 16 -S 0"
	echo "$EXEC"
	$EXEC || exit
fi



if true; then

	make clean

	SCONS="scons  --unit-test=test_antialiasing_patterns --plane-spectral-space=enable --mode=release --plane-spectral-dealiasing=enable"
	echo "$SCONS"
	$SCONS

	# test spectral derivatives
	EXEC="./build/test_antialiasing_patterns_*_release -N 16 -S 1"
	echo "$EXEC"
	$EXEC || exit

	# test cartesian derivatives
	EXEC="./build/test_antialiasing_patterns_*_release -N 16 -S 0"
	echo "$EXEC"
	$EXEC || exit

fi


echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
