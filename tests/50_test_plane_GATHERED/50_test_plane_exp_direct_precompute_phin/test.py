#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.exec_program import *
from mule.JobGeneration import *
from mule.SWEETRuntimeParametersScenarios import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_plane_exp_direct_precompute_phin"

earth = EarthMKSDimensions()
jg.runtime.verbosity = 6
jg.runtime.benchmark_name = "unstablejet"
jg = DisableGUI(jg)
jg.runtime.rexi_method = 'direct'
jg = RuntimeSWEPlaneEarthParam(jg)
jg.runtime.viscosity = 0.0
jg.runtime.max_simulation_time = 86400.
jg.runtime.timestepping_order = 2
jg.runtime.timestepping_order2 = 2
jg.runtime.space_res_physical = -1
jg.runtime.space_res_spectral = 512


jg.runtime.timestepping_method = "l_rexi_na_sl_nd_etdrk";
jg.runtime.timestep_size = 3600.;

jg.gen_jobscript_directory();


##exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
exitcode = os.system('job*/run.sh');
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)#! /usr/bin/env python3
