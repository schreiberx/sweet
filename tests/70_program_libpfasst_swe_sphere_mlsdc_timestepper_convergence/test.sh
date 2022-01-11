#! /bin/bash

cd "$(dirname $0)"

mule.benchmark.cleanup_all || exit 1

#./benchmark_create_job_scripts.py $TIMESTEPPING_GROUP || exit 1

#mule.benchmark.jobs_run_directly || exit 1

#./postprocessing_pickle.py || exit 1

#./postprocessing_convergence_test.py || exit 1

#mule.benchmark.cleanup_all || exit 1
