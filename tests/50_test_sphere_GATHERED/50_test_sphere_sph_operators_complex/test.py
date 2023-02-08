#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_sphere_sph_operators_complex"
jg.compile.plane_spectral_space="disable"
jg.compile.sphere_spectral_space="enable"
jg.compile.mode = 'debug'

params_runtime_mode_res = [64, 128, 256, 512, 1024, 2048]
jg.runtime.verbosity = 5

for jg.runtime.space_res_spectral in params_runtime_mode_res:
    jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
