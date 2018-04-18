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


# rotation speeds
# 0: no rotation
# 20=simtime: one rotation
# 1: fast rotaitons
for r in 0 1 20; do

	for o in 1 2; do
		for vu in 0.1 0.2 0.0; do
			EXEC="./build/test_plane_advection_*_debug -M 64 --dt=0.1 -X 2 -Y 2 --benchmark=gaussian_bump_advection --timestepping-method=na_sl --timestepping-order=$o  --advection-velocity=$vu,0.2,$r -t 20"
			echo "$EXEC"
			$EXEC || exit
		done
	done

	for vu in 0.1 0.2 0.0; do
		for o in 1 2 4; do
			EXEC="./build/test_plane_advection_*_debug -M 64 --dt=0.1 -X 2 -Y 2 --benchmark=gaussian_bump_advection --timestepping-method=na_erk --timestepping-order=$o  --advection-velocity=$vu,0.2,$r -t 20"
			echo "$EXEC"
			$EXEC || exit
		done
	done

done


echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
