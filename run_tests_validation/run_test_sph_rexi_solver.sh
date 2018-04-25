#! /bin/bash


echo "***********************************************"
echo "Running tests for SPH solver complex"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd ../

make clean
SCONS="scons --threading=omp --unit-test=test_sph_rexi_solver --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit 1

rm gen_*
#export SWEET_HACK_123_BLARG=1


for f in 0.1 10; do
#	for g in 0.1 10 100; do
#		for h in 0.1 10 100; do
	for g in 0.1 100; do
		for h in 0.1 100; do
			for r in 0.1 100; do
			    #EXEC="./build/test_sph_rexi_solver*_debug -f $f -g $g -H $h -a $r -M 256 --nonlinear 0 --use-robert-functions 1 --rexi-ext-modes 4 --rexi-m 2"
			    EXEC="./build/test_sph_rexi_solver*_debug -f $f -g $g -H $h -a $r -M 256 --use-robert-functions 1 --rexi-ext-modes 4 --rexi-m 2"
				echo "$EXEC"
				$EXEC || exit
			done
		done
	done
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
