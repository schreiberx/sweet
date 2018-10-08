#! /bin/bash


echo "***********************************************"
echo "Running geostrophic balance benchmark tests"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd "./run_test_swe_sphere_geostrophic_balance" || exit


./cleanup_all || exit 1

./benchmark_create_job_scripts || exit 1

./platform_jobs_run_directly || exit 1

./postprocessing || exit 1

./cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

