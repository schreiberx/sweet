#! /bin/bash


echo "***********************************************"
echo "Running convergence tests"
echo "***********************************************"

cd "run_test_swe_plane_timestepper_convergence"


./cleanup_all || exit 1

# encoding: group | tsm | order1 | order2 | rexi_direct_solution

# 2nd order nonlinear
./benchmark_create_job_scripts.py ln2 ln_erk 2 2 0 || exit 1
./benchmark_create_job_scripts.py ln2 l_cn_n_erk 2 2 0 || exit 1
./benchmark_create_job_scripts.py ln2 l_erk_n_erk 2 2 0 || exit 1
./benchmark_create_job_scripts.py ln2 l_rexi_n_erk 2 2 1 || exit 1
./benchmark_create_job_scripts.py ln2 l_rexi_n_etdrk 2 2 1 || exit 1


./platform_jobs_run_directly || exit 1

./postprocessing || exit 1

./cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

