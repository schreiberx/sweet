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

orders = {};
orders["l_irk_n_erk"] = 2;
orders["l_irn_na_sl_nd_settls"] = 2;

tsm_fine = "l_irk_n_erk"

#Create main compile/run options
jg = JobGeneration()


jg.runtime.output_file_mode = "csv"

# Verbosity mode
jg.runtime.verbosity = 3

jg.compile.fortran_source = "enable"
jg.compile.lapack = "enable"
jg.compile.sphere_spectral_space = "enable";
jg.compile.sphere_spectral_dealiasing = "enable";

#
# Benchmark ID
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "gaussian_bump"

#
# Compute error or difference to initial data
#
jg.runtime.compute_errors = 0

jg = DisableGUI(jg)

#
# REXI method
###jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1


jg.runtime.viscosity = 0.0

# Deactivate threading
jg.compile.threading = "omp"

#
# Time, Mode and Physical resolution
#
jg.runtime.max_simulation_time = 720.
jg.runtime.output_timestep_size = 180.
timestep_size_reference = 18.
timestep_size_fine = 36.
jg.runtime.timestep_size = timestep_size_fine
jg.runtime.timestepping_method = tsm_fine
jg.runtime.timestepping_order = orders[tsm_fine]
jg.runtime.timestepping_order2 = orders[tsm_fine]
jg.runtime.space_res_spectral = 32


## Reference job
jg.reference_job = True

jg.compile.program = "programs/pde_sweSphere"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "disable"

jg.compile.xbraid = "none";
jg.runtime.xbraid_enabled = 0;

jg.compile.parareal = "none";
jg.runtime.parareal_enabled = 0;

jg.runtime.timestepping_method = tsm_fine
jg.runtime.timestepping_order = orders[tsm_fine]
jg.runtime.timestepping_order2 = orders[tsm_fine]

jg.gen_jobscript_directory();


## MGRIT jobs

jg.reference_job = False
jg.reference_job_unique_id = jg.job_unique_id
ref_path = jg.p_job_dirpath
print("REF PATH ", ref_path)

jg.compile.program = "programs/xbraid_pde_sweSphere"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.compile.xbraid_sphere = 'enable'

jg.compile.xbraid = "mpi";
jg.runtime.xbraid_enabled = 1;
jg.runtime.xbraid_max_levels = 3
jg.runtime.xbraid_skip = 1
jg.runtime.xbraid_min_coarse = 2
jg.runtime.xbraid_nrelax = 1
jg.runtime.xbraid_nrelax0 = -1
jg.runtime.xbraid_tol = 1e-14
jg.runtime.xbraid_tnorm = 2
jg.runtime.xbraid_cfactor = 2
jg.runtime.xbraid_cfactor0 = -1
jg.runtime.xbraid_max_iter = 100
jg.runtime.xbraid_fmg = 0
jg.runtime.xbraid_res = 0
jg.runtime.xbraid_storage = 0
jg.runtime.xbraid_print_level = 2
jg.runtime.xbraid_access_level = 1
jg.runtime.xbraid_run_wrapper_tests = 0
jg.runtime.xbraid_fullrnorm = 2
jg.runtime.xbraid_use_seq_soln = 0
jg.runtime.xbraid_use_rand = 1
jg.runtime.xbraid_timestepping_method = tsm_fine
jg.runtime.xbraid_timestepping_order = orders[tsm_fine]
jg.runtime.xbraid_timestepping_order2 = orders[tsm_fine]
jg.runtime.xbraid_verbosity = 0;
jg.runtime.xbraid_load_ref_csv_files = 0;
jg.runtime.xbraid_path_ref_csv_files = "";
jg.runtime.xbraid_load_fine_csv_files = 0;
jg.runtime.xbraid_path_fine_csv_files = "";
jg.runtime.xbraid_store_iterations = 0;
jg.runtime.xbraid_spatial_coarsening = 0;
jg.runtime.xbraid_pt = 1;

jg.runtime.xbraid_store_iterations = 1;

nb_pts = [1]
##nb_pts = [1, 2, 4];

jg.runtime.xbraid_print_level = 3

for nb_pt in nb_pts:
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

