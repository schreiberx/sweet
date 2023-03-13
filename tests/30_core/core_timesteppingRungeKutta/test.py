#! /usr/bin/env python3

import sys
import os
import numpy as np

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()

jg.compile.program = "tests/core_timesteppingRungeKutta"
jg.runtime.verbosity = 5
jg.runtime.space_res_physical = [8,8]

for jg.runtime.max_simulation_time in [1, 2, 5, 7, np.pi]:
    for jg.runtime.timestep_size in [0.01, 0.05, 0.1, 0.123, np.pi*0.1]:
        jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
