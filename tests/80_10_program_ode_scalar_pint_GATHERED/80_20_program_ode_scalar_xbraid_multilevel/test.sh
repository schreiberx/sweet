#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> multilevel tests + n processors in time -> XBraid solution within tol w.r.t. serial solution
###############

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1

./benchmarks_create.py || exit 1

## force ref job to run first!
mule.benchmark.jobs_run_directly job_benchref* || exit 1

## run all jobs
mule.benchmark.jobs_run_directly job_bench_* || exit 1

mule.postprocessing.pickle.alljobs.scalar_data_norms_physical_space_pint

./postprocessing_check_errors.py

mule.benchmark.cleanup_job_dirs || exit 1
