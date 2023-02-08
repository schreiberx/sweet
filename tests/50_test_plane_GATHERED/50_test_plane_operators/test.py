#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_plane_operators"

jg.compile.plane_spectral_space="enable"

#params_domain_size_scales = [0.01, 1.0, 10000*1000]
params_domain_size_scales = [0.1, 10000*1000]
#params_domain_size_scales = [1.0]
#params_domain_size_scales = [10000*1000]

#params_compile_mode = ['release', 'debug']
params_compile_mode = ['debug']
params_compile_plane_spectral_dealiasing = ['enable', 'disable']

params_runtime_spectral_derivs = [0, 1]

params_runtime_phys_res_x = [16]
params_runtime_phys_res_y = [16]

for (res_x, res_y) in product(params_runtime_phys_res_x, params_runtime_phys_res_y):
    jg.runtime.space_res_physical = (res_x, res_y)

    # Try out different variants of domain size
    for jg.runtime.plane_domain_size in product(params_domain_size_scales, params_domain_size_scales):

    	for (
    		jg.compile.mode,
    		jg.compile.plane_spectral_dealiasing,
    		jg.runtime.space_use_spectral_basis_diffs,
    	) in product(
    		params_compile_mode,
    		params_compile_plane_spectral_dealiasing,
    		params_runtime_spectral_derivs
    	):
    		jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
