#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from SWEET import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = SWEETJobGeneration()
jg.compile.unit_test="test_sphere_sph_solver_real_and_complex"
jg.compile.plane_spectral_space="disable"
jg.compile.sphere_spectral_space="enable"
#jg.compile.mode = "debug"

jg.runtime.r = 1
jg.runtime.f = 1

params_runtime_mode_res = [64, 128, 256, 512, 1024, 2048]
#params_runtime_mode_res = [64, 128]

params_runtime_r = [1, 1e3, 1e6]
params_runtime_f = [1, 1e-3, 1e-6]

jg.runtime.verbosity = 5

for (
		jg.runtime.mode_res,
		jg.runtime.r,
		jg.runtime.f,
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
