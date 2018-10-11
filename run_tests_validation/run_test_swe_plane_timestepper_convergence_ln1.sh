#! /bin/bash


echo "***********************************************"
echo "Running convergence tests"
echo "***********************************************"

cd "run_test_swe_plane_timestepper_convergence"


./cleanup_all || exit 1

# encoding: group | tsm | order1 | order2 | rexi_direct_solution

# 1st order nonlinear
./benchmark_create_job_scripts.py ln1 ln_erk 1 1 0 || exit 1
./benchmark_create_job_scripts.py ln1 l_erk_n_erk 1 1 0 || exit 1
./benchmark_create_job_scripts.py ln1 l_irk_n_erk 1 1 0 || exit 1
./benchmark_create_job_scripts.py ln1 l_rexi_n_erk 1 1 1 || exit 1
./benchmark_create_job_scripts.py ln1 lg_rexi_lc_n_erk 1 1 1 || exit 1
./benchmark_create_job_scripts.py ln1 l_rexi_n_etdrk 1 1 1 || exit 1
./benchmark_create_job_scripts.py ln1 lg_rexi_lc_n_etdrk 1 1 1 || exit 1

./platform_jobs_run_directly || exit 1

./postprocessing || exit 1

./cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

