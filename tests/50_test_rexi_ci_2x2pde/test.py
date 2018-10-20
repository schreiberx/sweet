#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from SWEET import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = SWEETJobGeneration()

jg.compile.unit_test="test_rexi_ci_2x2pde"
jg.compile.quadmath="enable"
jg.runtime.verbosity=5

jg.runtime.rexi_method="ci"
jg.runtime.rexi_ci_primitive="circle"
jg.runtime.rexi_ci_n="64"

params_runtime_rexi_ci = [
			# n, real, imag
			[64, 10, 15],
			[96, 10, 25],
			[128, 10, 30],
			[256, 15, 50],
			[64, 7, 10],
]

jg.runtime.rexi_ci_sx = 25
jg.runtime.rexi_ci_sy = 25
jg.runtime.rexi_ci_n = 64

jg.gen_jobscript_directory()

jg.runtime.rexi_ci_sx = None
jg.runtime.rexi_ci_sy = None
jg.runtime.rexi_ci_n = None

for i in params_runtime_rexi_ci:
	jg.runtime.rexi_ci_n = i[0]
	jg.runtime.rexi_ci_max_real = i[1]
	jg.runtime.rexi_ci_max_imag = i[2]

	jg.gen_jobscript_directory()


exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
	sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
