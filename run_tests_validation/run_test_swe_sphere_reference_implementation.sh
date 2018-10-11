#! /bin/bash


echo "***********************************************"
echo "Running comparisons of SWEET SWE on the sphere implementation with Python reference implementation
echo "***********************************************"


cd "./run_test_swe_sphere_reference_implementation" || exit

echo "TODO"
./cleanup_all || exit 1

./benchmark_create_job_scripts || exit 1

./platform_jobs_run_directly || exit 1

./postprocessing || exit 1

./cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

