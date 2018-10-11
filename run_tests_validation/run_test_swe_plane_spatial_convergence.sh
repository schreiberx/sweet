#! /bin/bash


echo "***********************************************"
echo "Running convergence tests"
echo "***********************************************"

cd "run_test_swe_plane_spatial_convergence"

./cleanup_all || exit 1

./benchmark_create_jobs || exit 1

./platform_jobs_run_directly || exit 1

./postprocessing.py || exit 1

./cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

