#! /bin/bash


echo "***********************************************"
echo "Running comparisons of SWEET SWE on the sphere implementation with Python reference implementation"
echo "***********************************************"

cd "${SWEET_ROOT}"


mule.benchmark.cleanup_all || exit 1

./benchmark_create_job_scripts || exit 1

./benchmark_gen_swe_reference_solution.py || exit 1

mule.benchmark.jobs_run_directly || exit 1

./postprocessing || exit 1

mule.benchmark.cleanup_all || exit 1

echo "***********************************************"
echo "***************** FIN *************************"
echo "***********************************************"

