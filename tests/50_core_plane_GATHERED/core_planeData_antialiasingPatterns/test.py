#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from mule.utils import exec_program
from mule.InfoError import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="tests/core_planeData_antialiasingPatterns"
jg.compile.plane_spectral_space="enable"
jg.compile.plane_spectral_dealiasing="enable"
jg.compile.mode="release"
for jg.runtime.space_use_spectral_basis_diffs in [0, 1]:
    for nx in [16, 36]:
    	for ny in [16, 36]:
    		jg.runtime.space_res_physical = (nx, ny)
    		jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)

if exitcode != 0:
    sys.exit(exitcode)

exec_program('mule.benchmark.cleanup_all', catch_output=False)
