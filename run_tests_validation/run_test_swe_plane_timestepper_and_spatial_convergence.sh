#! /bin/bash


echo "***********************************************"
echo "Running convergence tests"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close


BASEDIR=`pwd`

cd "$BASEDIR/run_test_swe_timestepper_and_spatial_convergence"
./compile.sh || exit 1

./cleanup.sh


# 1st order linear
#./jobs_create.py l1 l_irk 1 0
#./jobs_create.py l1 l_erk 1 0

# 2nd order linear
#./jobs_create.py l2 l_cn 2 0
#./jobs_create.py l2 l_erk 2 0

# 1st order nonlinear
#./jobs_create.py ln1 ln_erk 1 1
#./jobs_create.py ln1 l_erk_n_erk 1 1
#./jobs_create.py ln1 l_irk_n_erk 1 1

# 2nd order nonlinear
#./jobs_create.py ln2 ln_erk 2 2
#./jobs_create.py ln2 l_cn_n_erk 2 2
#./jobs_create.py ln2 l_erk_n_erk 2 2

# 2nd order nonlinear, SL-REXI related
./jobs_create.py ln2space
#./jobs_create.py ln2 ln_erk 2 2

./jobs_run.sh || exit 1

echo "***********************************************"
echo " POSTPROCESSING "
echo "***********************************************"
./postprocessing_analytical.py || exit 1


#./cleanup.sh

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

