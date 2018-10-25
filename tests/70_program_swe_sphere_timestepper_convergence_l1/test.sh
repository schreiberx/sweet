#! /bin/bash

cd "$(dirname $0)"


TIMESTEPPING_GROUP="l1"

COMMON="../70_program_swe_sphere_timestepper_convergence_common/"


mule.benchmark.cleanup_all || exit 1

$COMMON/benchmark_create_job_scripts.py $TIMESTEPPING_GROUP || exit 1

mule.benchmark.jobs_run_directly || exit 1

$COMMON/postprocessing_pickle.py || exit 1

$COMMON/postprocessing_convergence_test.py || exit 1

mule.benchmark.cleanup_all || exit 1
