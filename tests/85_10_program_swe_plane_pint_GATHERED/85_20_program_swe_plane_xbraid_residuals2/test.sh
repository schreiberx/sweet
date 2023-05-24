#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> print_level = 3 -> check residual norm at each C-point: should be zero for two C-points per iteration
###############

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1

./benchmarks_create.py || exit 1

mule.benchmark.jobs_run_directly || exit 1

./check_residual.py C-point 1e-16

mule.benchmark.cleanup_job_dirs || exit 1
