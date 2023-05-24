#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> max_levels = 1 -> solution should be equal to the serial one
## --> max_levels = 1 and + processors in time -> idem
###############

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1

./benchmarks_create.py || exit 1

mule.benchmark.jobs_run_directly || exit 1

mule.postprocessing.pickle.alljobs.plane_data_norms_physical_space_pint

./postprocessing_check_errors.py

mule.benchmark.cleanup_job_dirs || exit 1
