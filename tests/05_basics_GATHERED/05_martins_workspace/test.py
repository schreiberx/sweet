#! /usr/bin/env python3

import sys
import os
import numpy as np

#Dummy email for tests on the linux cluster platform
os.environ['MULE_USER_EMAIL']='dummy_email@linuxcluster.com'

from mule.JobMule import *
from mule.utils import exec_program
from mule.InfoError import *


jg = JobGeneration()
jg.compile.mode = "debug"
jg.compile.program = 'programs/martin_workspace'

jg.runtime.user_defined_parameters['a-dict'] = {'id': '', 'value': "12", 'option': '--a='}
jg.runtime.user_defined_parameters['g-dict'] = {'id': '', 'value': "1", 'option': '-g='}
jg.runtime.user_defined_parameters['pde-viscosity'] = {'id': '', 'value': "12", 'option': '--pde-viscosity='}

jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)

