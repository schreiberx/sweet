#! /bin/bash


echo "***********************************************"
echo "Running convergence tests"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


BASEDIR=`pwd`

cd "$BASEDIR/run_test_swe_sphere_timestepper_convergence"


./cleanup.sh

# encoding: group | tsm | order1 | order2 | rexi_direct_solution

#if true; then
if false; then
	echo "Generating test scripts"

	#./jobs_create.py ln2 ln_erk 2 2 0 || exit 1
	./jobs_create.py ln2 l_rexi_n_etdrk 2 2 0 || exit 1
	#./jobs_create.py ln2 l_erk_n_erk 2 2 0 || exit 1

else
	# 2nd order nonlinear
	./jobs_create.py ln2 ln_erk 2 2 0 || exit 1
	./jobs_create.py ln2 l_rexi_n_etdrk 2 2 0 || exit 1
	./jobs_create.py ln2 lg_rexi_lc_n_etdrk 2 2 0 || exit 1
	#./jobs_create.py ln2 l_cn_n_erk 2 2 0 || exit 1
	./jobs_create.py ln2 l_erk_n_erk 2 2 0 || exit 1
	./jobs_create.py ln2 l_cn_n_erk 2 2 0 || exit 1
	#./jobs_create.py ln2 lg_cn_lc_n_erk 2 2 0 || exit 1
	#./jobs_create.py ln2 l_rexi_n_erk 2 2 1 || exit 1
fi


#./compile.sh || exit 1

./jobs_run.sh || exit 1



echo "***********************************************"
echo " POSTPROCESSING "
echo "***********************************************"
./postprocessing.py || exit 1


./cleanup.sh

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

