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

#for a in 1.5708 0 1.4 -0.7; do
for a in 0 1.4 -0.7; do
	EXEC="./build/test_plane_advection* -M 64 --dt=$((60*60*1)) -t $((60*60*24*12)) --benchmark=adv_gauss_bump --timestepping-order=2 --timestepping-method=na_sl --advection-rotation-angle=$a"
	echo "$EXEC"
	$EXEC || exit
done


for a in 1.5708 0 1.4 -0.7; do
	EXEC="./build/test_plane_advection* -M 64 --dt=$((60*60)) -t $((60*60*24*12)) --benchmark=adv_gauss_bump --timestepping-order=2 --timestepping-method=na_erk --advection-rotation-angle=$a"
	echo "$EXEC"
	#$EXEC || exit
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
