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

tsm_fine = "dummy";
tsm_coarse = "dummy";

#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
jg.compile.program = "parareal_ode"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"

jg.runtime.output_file_mode = "csv"

# Verbosity mode
jg.runtime.verbosity = 3

jg.compile.sphere_spectral_space = "enable";
jg.compile.sphere_spectral_dealiasing = "enable";

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "unstablejet"

#
# Compute error or difference to initial data
#
jg.runtime.compute_error = 1

# Enable/Disbale GUI
#jg = EnableGUI(jg)
jg = DisableGUI(jg)

#
# REXI method
jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1

# Parameters for SL-REXI paper
#-----------------------------
jg = RuntimeSWEPlaneEarthParam(jg)
#jg = RuntimeSWENonDimParam(jg)

jg.runtime.viscosity = 0.0

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
jg.runtime.space_res_spectral = 32
cfactors = [2, 4, 8];
nbs_levels = [2, 4];
nb_pts = [1];


## Reference job
jg.reference_job = True

jg.compile.xbraid = "none";
jg.runtime.xbraid_enabled = 0;

jg.compile.parareal = "none";
jg.runtime.parareal_enabled = 0;
jg.gen_jobscript_directory();


## MGRIT jobs

jg.reference_job = False
jg.reference_job_unique_id = jg.job_unique_id

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
jg.runtime.xbraid_timestepping_method = "dummy"
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

for nb_pt in range(1,4):
    jg.runtime.xbraid_pt = nb_pt;
    if nb_pt > 1:
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


