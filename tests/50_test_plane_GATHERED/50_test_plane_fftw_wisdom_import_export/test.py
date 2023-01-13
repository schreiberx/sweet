#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_plane_fftw_wisdom_import_export"
jg.compile.plane_spectral_space="enable"
jg.compile.plane_spectral_dealiasing="enable"
jg.compile.mode="debug"

jg.unique_id_filter = ['runtime.benchmark', 'runtime.time']

jg.runtime.space_res_physical = (128, 128)



if True:
    # No parallelization
    jg.compile.threading = 'off'

    pspace = JobParallelizationDimOptions('space')
    pspace.num_cores_per_rank = 1
    pspace.num_threads_per_rank = 1
    pspace.num_ranks = 1
    jg.setup_parallelization(pspace)

    # Create plans
    jg.runtime.reuse_plans = "save"
    jg.gen_jobscript_directory()

    # Reuse plans
    jg.runtime.reuse_plans = "require_load"
    jg.gen_jobscript_directory()




if True:
    # Parallelization
    jg.compile.threading = 'omp'
    for i in range(1, jg.platform_resources.num_cores_per_socket+1):
    	pspace = JobParallelizationDimOptions('space')
    	pspace.num_cores_per_rank = i
    	pspace.num_ranks = 1
    	jg.setup_parallelization(pspace)

    	# Create plans
    	jg.runtime.reuse_plans = "save"
    	jg.gen_jobscript_directory()

    	# Reuse plans
    	jg.runtime.reuse_plans = "require_load"
    	jg.gen_jobscript_directory()



import glob

print("Running jobs...")

for i in glob.glob("job_bench_*_planssave_*"):
    print("Executing: "+i)
    exitcode = exec_program(['mule.benchmark.jobs_run_directly', i], catch_output=False)
    if exitcode != 0:
        sys.exit(exitcode)

for i in glob.glob("job_bench_*_plansrequire_load_*"):
    print("Executing: "+i)
    exitcode = exec_program(['mule.benchmark.jobs_run_directly', i], catch_output=False)
    if exitcode != 0:
        sys.exit(exitcode)

exec_program('mule.benchmark.cleanup_all', catch_output=False)
