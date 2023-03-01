#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> max_levels = 2, tol = 0., max_iter = 3, use_seqsoln, n processors in time -> residual should be zero at each iteration
###############

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1

./benchmarks_create.py || exit 1

mule.benchmark.jobs_run_directly || exit 1

./check_residual.py iteration 1e-16

mule.benchmark.cleanup_job_dirs || exit 1
