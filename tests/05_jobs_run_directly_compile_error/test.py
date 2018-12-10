#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from itertools import product

from mule_local.JobMule import *
from mule.exec_program import *

ie = InfoError("test")

exec_program('mule.benchmark.cleanup_all', catch_output=False)


"""
Generate a job which fails.
This assures that we are able to catch failing jobs in case of malfunctioning scripts
"""

jg = JobGeneration()
jg.compile.program="I_DONT_EXIST"
jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode == 0:
	ie.error("The job should return with exit code != 0")
	ie.error("Benchmark's return value != 0")
	sys.exit(1)

ie.success_hline()
ie.success_hline()
ie.success_hline()
ie.success("")
ie.success("The benchmark job script was intended to fail!")
ie.success("Benchmarks successfully finished")
ie.success("")
ie.success_hline()
ie.success_hline()
ie.success_hline()

exec_program('mule.benchmark.cleanup_all', catch_output=False)
sys.exit(0)
