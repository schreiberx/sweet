#! /bin/bash


echo "***********************************************"
echo "Running geostrophic balance benchmark tests"
echo "***********************************************"

# set close affinity of threads
export OMP_PROC_BIND=close

cd "./run_test_swe_sphere_geostrophic_balance" || exit


./cleanup_all

./jobs_create_scripts

./platform_jobs_run || exit 1

./postprocessing || exit 1

./cleanup_all

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

