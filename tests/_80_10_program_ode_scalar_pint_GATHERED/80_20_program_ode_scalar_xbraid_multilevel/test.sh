#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> multilevel tests + n processors in time -> XBraid solution within tol w.r.t. serial solution
###############

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1

./benchmarks_create.py || exit 1

## force ref job to be run first!
mkdir tmp
mv job_bench_* tmp/.
mule.benchmark.jobs_run_directly || exit 1

## run all jobs
mv tmp/* .
rm -r tmp
mule.benchmark.jobs_run_directly || exit 1

mule.postprocessing.pickle.alljobs.scalar_data_norms_physical_space_pint

./postprocessing_check_errors.py

mule.benchmark.cleanup_job_dirs || exit 1
