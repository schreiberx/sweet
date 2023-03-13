#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="tests/core_memBlockAlloc"
jg.compile.mode="debug"

jg.unique_id_filter = ['runtime.benchmark', 'runtime.time']


jg.compile.threading = "omp"

jg.gen_jobscript_directory()


exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

exec_program('mule.benchmark.cleanup_all', catch_output=False)
