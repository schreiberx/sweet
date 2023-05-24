#! /usr/bin/env python3
#
#  Create series of job to be run with sweet
#
#  Pedro Peixoto <pedrosp@ime.usp.br>
#  modified from Martin Schreiber initial job_create
#
#-------------------------------------------------------

import os
import sys
import stat
import math
from glob import glob

#Classes containing sweet compile/run basic option
from mule.JobGeneration import *
from mule.SWEETRuntimeParametersScenarios import *

jg = JobGeneration()

#Get Earth parameters (if necessary)
#
# Run simulation on plane or sphere
#
#Basic plane options
jg.compile.program = "programs/xbraid_ode_Scalar"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.compile.xbraid_scalar = "enable"

jg.runtime.output_file_mode = "csv"

# Verbosity mode
jg.runtime.verbosity = 3

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "default"
jg.runtime.u0 = 1.
jg.runtime.param_a = 1.
jg.runtime.param_b = 2.

#
# Compute error or difference to initial data
#
jg.runtime.compute_errors = 1

# Enable/Disbale GUI
#jg = EnableGUI(jg)
jg = DisableGUI(jg)

# Deactivate threading
jg.compile.threading = 'off'


#
# Time, Mode and Physical resolution
#
jg.runtime.max_simulation_time = 1.
jg.runtime.output_timestep_size = 1.
jg.runtime.timestep_size = 1.

jg.compile.xbraid = "mpi"
jg.runtime.xbraid_enabled = 1
jg.runtime.xbraid_max_levels = 3
jg.runtime.xbraid_skip = 0
jg.runtime.xbraid_min_coarse = 2
jg.runtime.xbraid_nrelax = 1
jg.runtime.xbraid_nrelax0 = -1
jg.runtime.xbraid_tol = 1e-14
jg.runtime.xbraid_tnorm = 2
jg.runtime.xbraid_cfactor = 2
jg.runtime.xbraid_cfactor0 = -1
jg.runtime.xbraid_max_iter = 10
jg.runtime.xbraid_fmg = 0
jg.runtime.xbraid_res = 0
jg.runtime.xbraid_storage = 0
jg.runtime.xbraid_print_level = 2
jg.runtime.xbraid_access_level = 1
jg.runtime.xbraid_fullrnorm = 2
jg.runtime.xbraid_use_seq_soln = 0
jg.runtime.xbraid_use_rand = 1
jg.runtime.xbraid_timestepping_method = "ln_erk"
jg.runtime.xbraid_timestepping_order = "2"
jg.runtime.xbraid_timestepping_order2 = "2"
jg.runtime.xbraid_verbosity = 0
jg.runtime.xbraid_load_ref_csv_files = 0
jg.runtime.xbraid_path_ref_csv_files = ""
jg.runtime.xbraid_load_fine_csv_files = 0
jg.runtime.xbraid_path_fine_csv_files = ""
jg.runtime.xbraid_store_iterations = 0

jg.runtime.xbraid_run_wrapper_tests = 1

jg.gen_jobscript_directory()


