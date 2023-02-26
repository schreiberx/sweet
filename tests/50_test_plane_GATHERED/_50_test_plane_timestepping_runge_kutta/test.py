#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()

jg.compile.program = "tests/core_timesteppingRungeKutta"
jg.runtime.verbosity = 5
jg.runtime.max_simulation_time = 5
jg.runtime.max_simulation_time = 5
jg.runtime.timestep_size = 0.01

jg.runtime.space_res_physical = [8,8]

jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
