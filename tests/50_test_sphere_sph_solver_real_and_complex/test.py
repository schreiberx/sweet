#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule_local.JobMule import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_sphere_sph_solver_real_and_complex"

jg.compile.plane_spectral_space="disable"
jg.compile.sphere_spectral_space="enable"
jg.compile.mode = "release"

jg.runtime.sphere_radius = 1
jg.runtime.sphere_rotating_coriolis_omega = 1

unique_id_filter = []
unique_id_filter.append('compile')


jg.unique_id_filter = unique_id_filter


#params_runtime_mode_res = [64, 128, 256, 512, 1024, 2048]
params_runtime_mode_res = [64, 128, 256, 512, 1024]

params_runtime_r = [1, 1e3, 1e6]
params_runtime_f = [1, 1e-3, 1e-6]

jg.runtime.verbosity = 5

for (
    	jg.runtime.space_res_spectral,
    	jg.runtime.sphere_radius,
    	jg.runtime.sphere_rotating_coriolis_omega,
    ) in product(
    	params_runtime_mode_res,
    	params_runtime_r,
    	params_runtime_f,
    ):
    jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)

