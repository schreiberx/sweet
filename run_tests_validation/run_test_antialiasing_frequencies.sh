#! /bin/bash


echo "***********************************************"
echo "Running tests for antialiasing frequencies"
echo "***********************************************"

cd ../


for MODE in debug release; do
	make clean

	SCONS="scons  --unit-test=test_antialiasing_frequencies --plane-spectral-space=enable --mode=$MODE --plane-spectral-dealiasing=enable"
	echo "$SCONS"
	$SCONS

	for Nx in `seq 4 8 36`; do
		for Ny in `seq 4 8 36`; do
			EXEC="./build/test_antialiasing_frequencies_*_$MODE -N $Nx,$Ny"
			echo "$EXEC"
			$EXEC || exit
		done
	done

	EXEC="./build/test_antialiasing_frequencies_*_$MODE -N 64,64"
	echo "$EXEC"
	$EXEC || exit
done


echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
echo "***********************************************"
