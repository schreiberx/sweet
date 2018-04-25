#! /bin/bash


echo "***********************************************"
echo "Running tests for advection schemes on the sphere"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_sphere_advection --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=debug"
#SCONS="scons --threading=omp --unit-test=test_sphere_advection --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=release"
echo "$SCONS"
$SCONS || exit 1





for a in 1.5708 0 1.4 -0.7; do
	for r in 0 $((60*60*24*12)); do
		EXEC="./build/test_sphere_advection* -M 64 --dt=$((60*60*4)) -t $((60*60*24*12)) --benchmark=adv_gauss_bump --timestepping-order=2 --timestepping-method=na_sl --advection-rotation-angle=$a --advection-velocity=0,0,$r"
		echo "$EXEC"
		$EXEC || exit
		echo
		echo
		echo
	done
done

for a in 1.5708 0 1.4 -0.7; do
#for a in 0 1.4 -0.7; do
	EXEC="./build/test_sphere_advection* -M 64 --dt=$((60*60*6)) -t $((60*60*24*12)) --benchmark=adv_gauss_bump --timestepping-order=1 --timestepping-method=na_sl --advection-rotation-angle=$a"
	echo "$EXEC"
	$EXEC || exit
	echo
	echo
	echo
done


for a in 1.5708 0 1.4 -0.7; do
	EXEC="./build/test_sphere_advection* -M 64 --dt=$((60*60)) -t $((60*60*24*12)) --benchmark=adv_gauss_bump --timestepping-order=2 --timestepping-method=na_erk --advection-rotation-angle=$a"
	echo "$EXEC"
	$EXEC || exit
	echo
	echo
	echo
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
