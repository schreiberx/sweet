#! /bin/bash

cd "$(dirname $0)"


TIMESTEPPING_GROUP="ln1"


mule.benchmark.cleanup_all || exit 1

../50_test_plane_timestepper_convergence_common/benchmark_create_job_scripts.py $TIMESTEPPING_GROUP || exit 1

mule.benchmark.jobs_run_directly || exit 1

../50_test_plane_timestepper_convergence_common/postprocessing_pickle.py || exit 1

../50_test_plane_timestepper_convergence_common/postprocessing_convergence_test.py || exit 1

mule.benchmark.cleanup_all || exit 1
