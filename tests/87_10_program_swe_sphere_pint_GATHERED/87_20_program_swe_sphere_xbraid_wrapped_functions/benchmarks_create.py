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
tsm_coarse = "l_irk_n_erk"

#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
jg.compile.program = "programs/xbraid_pde_sweSphere"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.compile.xbraid_sphere = 'enable'
jg.compile.fortran_source = "enable"
jg.compile.lapack = "enable"

jg.runtime.output_file_mode = "csv"

# Verbosity mode
jg.runtime.verbosity = 3

jg.compile.sphere_spectral_space = "enable";
jg.compile.sphere_spectral_dealiasing = "enable";

#
# Benchmark ID
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "galewsky"

#
# Compute error or difference to initial data
#
jg.runtime.compute_errors = 1

jg = DisableGUI(jg)

#
# REXI method
###jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1

#jg = RuntimeSWEPlaneEarthParam(jg)
#jg = RuntimeSWENonDimParam(jg)

jg.runtime.viscosity = 0.0

# Deactivate threading
jg.compile.threading = "omp"

#
# Time, Mode and Physical resolution
#
jg.runtime.max_simulation_time = 100.
jg.runtime.output_timestep_size = 25.
timestep_size_reference = 2.5
timestep_size_fine = 5.
jg.runtime.timestep_size = timestep_size_fine
jg.runtime.timestepping_method = tsm_fine
jg.runtime.timestepping_order = orders[tsm_fine]
jg.runtime.timestepping_order2 = orders[tsm_fine]
jg.runtime.space_res_spectral = 32

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

jg.runtime.xbraid_run_wrapper_tests = 1

jg.gen_jobscript_directory();
