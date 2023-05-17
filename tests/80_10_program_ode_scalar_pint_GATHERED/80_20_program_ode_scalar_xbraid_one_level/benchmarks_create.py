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
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *

tsm_fine = "ln_erk";
tsm_coarse = "ln_erk";

#Create main compile/run options
jg = JobGeneration()


jg.runtime.output_file_mode = "csv"

# Verbosity mode
jg.runtime.verbosity = 3


#
# Benchmark ID
#
jg.runtime.benchmark_name = "default"
jg.runtime.u0 = 1.
jg.runtime.param_a = 1.
jg.runtime.param_b = 2.

#
# Compute error or difference to initial data
#
jg.runtime.compute_errors = 1

# Enable/Disbale GUI
jg = DisableGUI(jg)


# Deactivate threading
jg.compile.threading = 'off'

#
# Time, Mode and Physical resolution
#
jg.runtime.max_simulation_time = 1.
jg.runtime.output_timestep_size = .1

timestep_size_reference = 0.001
timestep_size_fine = 0.005
jg.runtime.timestep_size = timestep_size_fine
cfactors = [2, 4, 8];
nbs_levels = [2, 4];
nb_pts = [1];


## Reference job
jg.reference_job = True

jg.compile.program = "programs/ode_Scalar"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "disable"
jg.compile.xbraid_scalar = 'disable'

jg.compile.xbraid = "none";
jg.runtime.xbraid_enabled = 0;

jg.compile.parareal = "none";
jg.runtime.parareal_enabled = 0;

jg.runtime.timestepping_method = "ln_erk"
jg.runtime.timestepping_order = 2
jg.runtime.timestepping_order2 = 2

jg.gen_jobscript_directory();


## MGRIT jobs

jg.reference_job = False
jg.reference_job_unique_id = jg.job_unique_id

jg.compile.program = "programs/xbraid_ode_Scalar"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.compile.xbraid_scalar = 'enable'

jg.compile.xbraid = "mpi";
jg.runtime.xbraid_enabled = 1;
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
jg.runtime.xbraid_run_wrapper_tests = 0
jg.runtime.xbraid_fullrnorm = 2
jg.runtime.xbraid_use_seq_soln = 0
jg.runtime.xbraid_use_rand = 1
jg.runtime.xbraid_timestepping_method = "ln_erk"
jg.runtime.xbraid_timestepping_order = "2"
jg.runtime.xbraid_timestepping_order2 = "2"
jg.runtime.xbraid_verbosity = 0;
jg.runtime.xbraid_load_ref_csv_files = 0;
jg.runtime.xbraid_path_ref_csv_files = "";
jg.runtime.xbraid_load_fine_csv_files = 0;
jg.runtime.xbraid_path_fine_csv_files = "";
jg.runtime.xbraid_store_iterations = 0;

jg.runtime.xbraid_max_levels = 1
jg.runtime.xbraid_store_iterations = 1;

for nb_pt in range(1,5):
    jg.runtime.xbraid_pt = nb_pt;

    params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
    params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
    params_ptime_num_cores_per_rank = [1]

    # Update TIME parallelization
    ptime = JobParallelizationDimOptions('time')
    ptime.num_cores_per_rank = 1
    ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
    ptime.num_ranks = nb_pt

    pspace = JobParallelizationDimOptions('space')
    pspace.num_cores_per_rank = 1
    pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
    pspace.num_ranks = 1

    # Setup parallelization
    jg.setup_parallelization([pspace, ptime], override_insufficient_resources=True)

    jg.gen_jobscript_directory();


