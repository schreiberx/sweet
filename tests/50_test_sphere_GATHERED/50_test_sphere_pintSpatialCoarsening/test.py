#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.mode = "debug"
jg.compile.xbraid_sphere = "enable"
jg.compile.program="tests/core_sphere_pintSpatialCoarsening"
jg.compile.plane_spectral_space="disable"
jg.compile.sphere_spectral_space="enable"


jg.runtime.benchmark_name = "rossby_haurwitz_wave"
jg.runtime.space_res_spectral = 256;
jg.runtime.verbosity = 0

jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
