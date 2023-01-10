#! /bin/bash

cd "$(dirname $0)"


TIMESTEPPING_GROUPS="l1 l2 ln1 ln2"

COMMON="./"


mule.benchmark.cleanup_all || exit 1

for TIMESTEPPING_GROUP in $TIMESTEPPING_GROUPS; do
	$COMMON/benchmark_create_job_scripts.py $TIMESTEPPING_GROUP || exit 1
done

mule.benchmark.jobs_run_directly || exit 1

$COMMON/postprocessing_pickle.py || exit 1

$COMMON/postprocessing_convergence_test.py || exit 1

mule.benchmark.cleanup_all || exit 1

