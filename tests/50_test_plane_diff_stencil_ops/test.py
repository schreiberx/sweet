#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from SWEET import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = SWEETJobGeneration()
jg.compile.unit_test="test_plane_operators_diff_ops"

jg.compile.plane_spectral_space="enable"

params_domain_size_scales = [0.01, 10000*1000]

params_compile_mode = ['release', 'debug']
params_compile_plane_spectral_dealiasing = ['enable', 'disable']

#params_runtime_spectral_derivs = [0, 1]
params_runtime_spectral_derivs = [0]

params_runtime_phys_res_x = [128]
params_runtime_phys_res_y = [128]

params_runtime_mode_res_x = [-1]
params_runtime_mode_res_y = [-1]

for (phys_res_x, phys_res_y, mode_res_x, mode_res_y) in product(params_runtime_phys_res_x, params_runtime_phys_res_y, params_runtime_mode_res_x, params_runtime_mode_res_y):
	jg.runtime.phys_res = (phys_res_x, phys_res_y)
	jg.runtime.mode_res = (mode_res_x, mode_res_y)

	# Try out different variants of domain size
	for jg.runtime.domain_size in product(params_domain_size_scales, params_domain_size_scales):
#	if True:

		for (
			jg.compile.mode,
			jg.compile.plane_spectral_dealiasing,
			jg.runtime.spectralderiv,
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
