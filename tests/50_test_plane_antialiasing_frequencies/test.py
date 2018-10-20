#! /usr/bin/env python3

import sys

from SWEET import *
from mule.exec_program import *
from InfoError import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = SWEETJobGeneration()
jg.compile.unit_test="test_plane_antialiasing_frequencies"
jg.compile.plane_spectral_space="enable"
jg.compile.plane_spectral_dealiasing="enable"
jg.compile.mode="release"
for nx in [4, 8, 36]:
	for ny in [4, 8, 36]:
		jg.runtime.phys_res = (nx, ny)
		jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
	sys.exit(exitcode)

exec_program('mule.benchmark.cleanup_all', catch_output=False)
