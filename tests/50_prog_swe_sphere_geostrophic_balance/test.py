#! /usr/bin/env python3

import sys
import os

from SWEET import *
from mule.exec_program import *
from InfoError import *

os.chdir(os.path.dirname(sys.argv[0]))

exec_program('./benchmark_create_job_scripts', catch_output=False)

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
	sys.exit(exitcode)

exitcode = exec_program('./postprocessing', catch_output=False)
if exitcode != 0:
	sys.exit(exitcode)

exec_program('mule.benchmark.cleanup_all', catch_output=False)
