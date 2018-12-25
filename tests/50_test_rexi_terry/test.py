#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule_local.JobMule import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()

jg.compile.unit_test="test_rexi_terry"
jg.compile.quadmath="enable"
jg.runtime.verbosity=5

jg.runtime.rexi_method="terry"
jg.runtime.rexi_terry_m=256

jg.runtime.rexi_terry_reduce_to_half=0
jg.gen_jobscript_directory()

jg.runtime.rexi_terry_reduce_to_half=1
jg.gen_jobscript_directory()



exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
	sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
