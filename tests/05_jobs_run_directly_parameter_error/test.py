#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from SWEET import *
from itertools import product
from mule.exec_program import *

ie = InfoError("test")

exec_program('mule.benchmark.cleanup_all', catch_output=False)

"""
Generate a job which fails.
This assures that we are able to catch failing jobs in case of malfunctioning scripts
"""

jg = SWEETJobGeneration()
jg.compile.program="swe_plane"
jg.runtime.space_res_physical=128
jg.runtime.timestepping_method="ln_erk"
jg.runtime.timestepping_order=2
jg.runtime.timestep_size = 0.0001
jg.runtime.max_timesteps = 10
jg.runtime.benchmark_name="I_DONT_EXIST"	# Use benchmark which doesn't exist
jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode == 0:
	ie.error("Benchmark's return value != 0")
	sys.exit(1)

ie.success_hline()
ie.success("The benchmark job script was intended to fail!")
ie.success("Benchmarks successfully finished")
ie.success_hline()

exec_program('mule.benchmark.cleanup_all', catch_output=False)
sys.exit(0)
