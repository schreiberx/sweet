#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="tests/core_planeData_convertRealToComplex"

jg.compile.plane_spectral_space="enable"

params_compile_mode = ['release', 'debug']
params_compile_plane_spectral_dealiasing = ['enable', 'disable']

params_runtime_phys_res_x = [-1]
params_runtime_phys_res_y = [-1]

params_runtime_mode_res_x = [16]
params_runtime_mode_res_y = [16]

for (phys_res_x, phys_res_y, mode_res_x, mode_res_y) in product(params_runtime_phys_res_x, params_runtime_phys_res_y, params_runtime_mode_res_x, params_runtime_mode_res_y):
    jg.runtime.space_res_physical = (phys_res_x, phys_res_y)
    jg.runtime.space_res_physical = (mode_res_x, mode_res_y)

    # Try out different variants of domain size
#    for jg.runtime.plane_domain_size in product(params_domain_size_scales, params_domain_size_scales):
    if True:

    	for (
    		jg.compile.mode,
    		jg.compile.plane_spectral_dealiasing,
    	) in product(
    		params_compile_mode,
    		params_compile_plane_spectral_dealiasing,
    	):
    		jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
