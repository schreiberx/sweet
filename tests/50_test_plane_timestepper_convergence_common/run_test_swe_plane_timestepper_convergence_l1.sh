#! /bin/bash


echo "***********************************************"
echo "Running convergence tests"
echo "***********************************************"

cd "run_test_swe_plane_timestepper_convergence"


./cleanup_all || exit 1

# encoding: group | tsm | order1 | order2 | rexi_direct_solution
./benchmark_create_job_scripts.py l1 l_irk 1 0 0 || exit 1
./benchmark_create_job_scripts.py l1 l_erk 1 0 0 || exit 1

./platform_jobs_run_directly || exit 1

./postprocessing || exit 1

./cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

