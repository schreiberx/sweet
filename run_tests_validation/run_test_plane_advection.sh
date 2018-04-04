#! /bin/bash


echo "***********************************************"
echo "Running tests for advection schemes on the plane"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
#SCONS="scons --threading=omp --unit-test=plane_advection --gui=disable --plane-spectral-space=disable --plane-spectral-space=enable --mode=debug"
SCONS="scons --threading=omp --unit-test=test_plane_advection --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS


#for vu in 0.1 0.2 0.0; do
#	EXEC="./build/test_plane_advection_*_debug -M 64 --dt=0.1 -X 2 -Y 2 --benchmark=gaussian_bump --timestepping-method=na_erk --timestepping-order=4  --velocity-u=$vu --velocity-v=0.2 -t 20"
#	echo "$EXEC"
#	$EXEC || exit
#done

for vu in 0.1 0.2 0.0; do
	EXEC="./build/test_plane_advection_*_debug -M 64 --dt=0.1 -X 2 -Y 2 --benchmark=gaussian_bump --timestepping-method=na_sl --timestepping-order=4  --velocity-u=$vu --velocity-v=0.2 -t 20"
	echo "$EXEC"
	$EXEC || exit
done




echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
