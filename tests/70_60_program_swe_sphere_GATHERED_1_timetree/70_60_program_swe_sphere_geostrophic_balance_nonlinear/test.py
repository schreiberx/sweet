#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from mule.utils import exec_program
from mule.InfoError import *

exitcode = exec_program('./benchmark_create_job_scripts.py', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

exitcode = exec_program('./postprocessing.py', catch_output=False)
if exitcode != 0:
    print("FAILED")
    sys.exit(exitcode)

exec_program('mule.benchmark.cleanup_all', catch_output=False)
