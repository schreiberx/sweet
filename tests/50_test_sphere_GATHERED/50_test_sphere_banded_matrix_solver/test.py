#! /usr/bin/env python3

import sys
import os

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_banded_matrix_solver"
jg.compile.mode = 'debug'

jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
