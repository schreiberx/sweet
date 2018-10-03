#! /bin/bash


echo "***********************************************"
echo "Running tests for advection schemes on the plane"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
#SCONS="scons --threading=omp --unit-test=test_plane_advection --gui=disable --plane-spectral-space=enable --mode=debug"
SCONS="scons --threading=omp --unit-test=test_plane_advection --gui=disable --plane-spectral-space=enable --mode=release"
echo "$SCONS"
$SCONS || exit 1




#./build/test_plane_advection_* -M 64 --dt=0.01 -X 2 -Y 2 --benchmark-name=gaussian_bump_advection --timestepping-method=na_erk --timestepping-order=4  --advection-velocity=0.1,0.2,1 -t 20


#for r in 0 1 20; do
for r in 1 20; do
	for vu in 0.1 -0.2; do
		for vv in -0.1 0.2; do
			# order 1
			EXEC="./build/test_plane_advection_* -M 128 --dt=0.05 -X 2 -Y 2 --benchmark-name=gaussian_bump_advection --timestepping-method=na_sl --timestepping-order=1  --advection-velocity=$vu,$vv,$r -t 20"
			echo "$EXEC"
			$EXEC || exit
		done
	done
done


# rotation speeds
# 0: no rotation
# 20=simtime: one rotation
# 1: fast rotaitons
for r in 0 1 20; do
	for vu in 0.1 -0.2; do
		for vv in -0.1 0.2; do
			# order 2
			EXEC="./build/test_plane_advection_* -M 64 --dt=0.1 -X 2 -Y 2 --benchmark-name=gaussian_bump_advection --timestepping-method=na_sl --timestepping-order=2  --advection-velocity=$vu,$vv,$r -t 20"
			echo "$EXEC"
			$EXEC || exit
		done
	done
done


for r in 0 20; do
	for vu in 0.1 -0.2; do
		for vv in -0.1 0.2; do
			# order 4
			EXEC="./build/test_plane_advection_* -M 64 --dt=0.02 -X 2 -Y 2 --benchmark-name=gaussian_bump_advection --timestepping-method=na_erk --timestepping-order=4  --advection-velocity=$vu,$vv,$r -t 20"
			echo "$EXEC"
			$EXEC || exit

			# order 2
			EXEC="./build/test_plane_advection_* -M 64 --dt=0.01 -X 2 -Y 2 --benchmark-name=gaussian_bump_advection --timestepping-method=na_erk --timestepping-order=2  --advection-velocity=$vu,$vv,$r -t 20"
			echo "$EXEC"
			$EXEC || exit

			# order 1
#			EXEC="./build/test_plane_advection_* -M 64 --dt=0.002 -X 2 -Y 2 --benchmark-name=gaussian_bump_advection --timestepping-method=na_erk --timestepping-order=1  --advection-velocity=$vu,$vv,$r -t 20"
#			echo "$EXEC"
#			$EXEC || exit

		done
	done

done


echo "***********************************************"
echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"
echo "***********************************************"
