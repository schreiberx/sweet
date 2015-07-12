#! /bin/bash


echo "***********************************************"
echo "Running tests for validation of SWEET"
echo "***********************************************"



if true; then
	X=10000000

	echo
	echo "***********************************************"
	echo "TEST SPECTRAL OPS (debug) ALIASING CONTROL $X x $X"
	echo "***********************************************"
	make clean
	scons --spectral-dealiasing=enable --compile-program=test_spectral_ops --gui=disable --spectral-space=enable --mode=debug --spectral-dealiasing=disable
	./build/example_test_spectral_ops_spectral_gnu_debug -n 128 -m 64 -X $X -Y $X -S 0 || exit
	./build/example_test_spectral_ops_spectral_gnu_debug -n 128 -m 64 -X $X -Y $X -S 1 || exit


	echo
	echo "***********************************************"
	echo "TEST SPECTRAL OPS (debug) $X x $X"
	echo "***********************************************"
	make clean
	scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable --mode=debug --spectral-dealiasing=disable
	./build/example_test_spectral_ops_spectral_gnu_debug  -X $X -Y $X -S 0 || exit
	./build/example_test_spectral_ops_spectral_gnu_debug  -X $X -Y $X -S 1 || exit

	echo
	echo "***********************************************"
	echo "TEST SPECTRAL OPS (release) $X x $X"
	echo "***********************************************"
	make clean
	scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable --mode=release --spectral-dealiasing=disable
	./build/example_test_spectral_ops_spectral_gnu_release  -X $X -Y $X -S 0 || exit
	./build/example_test_spectral_ops_spectral_gnu_release  -X $X -Y $X -S 1 || exit


	echo
	echo "***********************************************"
	echo "TEST SPECTRAL OPS (release) $X x $X"
	echo "***********************************************"
	make clean
	scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable --mode=release --spectral-dealiasing=disable
	./build/example_test_spectral_ops_spectral_gnu_release -n 128 -m 64 -X $X -Y 500 -S 0 || exit
	./build/example_test_spectral_ops_spectral_gnu_release -n 128 -m 64 -X $X -Y 500 -S 1 || exit

	echo
	echo "***********************************************"
	echo "TEST SPECTRAL OPS (debug)"
	echo "***********************************************"
	make clean
	scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable --mode=debug --spectral-dealiasing=disable
	./build/example_test_spectral_ops_spectral_gnu_debug  -X 1 -Y 1 -S 0 || exit
	./build/example_test_spectral_ops_spectral_gnu_debug  -X 1 -Y 1 -S 1 || exit

	echo
	echo "***********************************************"
	echo "TEST SPECTRAL OPS (release)"
	echo "***********************************************"
	make clean
	scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable --mode=release --spectral-dealiasing=disable
	./build/example_test_spectral_ops_spectral_gnu_release  -X 1 -Y 1 -S 0 || exit
	./build/example_test_spectral_ops_spectral_gnu_release  -X 1 -Y 1 -S 1 || exit
fi
