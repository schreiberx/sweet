#! /bin/bash


#cd ../../ || exit 1

make clean
scons --program=libpfasst_swe_plane --libpfasst=enable --sweet-mpi=enable --pfasst-cpp=enable --libfft=enable || exit 1 

echo
echo "*********************************************************"
echo "* Running benchmarks                                    *"
echo "*********************************************************"
echo

EXEC="./build/libpfasst_swe_plane_* -m 8 -n 8 -M 8"
	echo "$EXEC"
	valgrind --leak-check=full $EXEC | grep "DIAGNOSTICS ANALYTICAL MAXABS H"
