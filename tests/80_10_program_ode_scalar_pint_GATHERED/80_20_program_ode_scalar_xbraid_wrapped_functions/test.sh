#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> Test user-defined wrapped functions
###############

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1

./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse 1 > tmp_job_benchmark_create_dummy.txt || exit 1

mule.benchmark.jobs_run_directly || exit 1

mule.benchmark.cleanup_job_dirs || exit 1
