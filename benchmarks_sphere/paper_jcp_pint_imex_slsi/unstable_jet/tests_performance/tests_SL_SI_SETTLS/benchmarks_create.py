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

#Create main compile/run options
jg = JobGeneration()

if os.environ['MULE_PLATFORM_ID'] == 'dahu_llvm':
    nb_cores_per_rank = 32
    nb_cores_per_rank_pt1 = 2
    nb_threads = 16
    nb_pts = [1, 2, 4, 8, 16]; ## number of parallel processors in time
else:
    nb_cores_per_rank = os.cpu_count()
    nb_cores_per_rank_pt1 = 1
    nb_threads = os.cpu_count()
    nb_pts = [1]; ## number of parallel processors in time


jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"


jg.compile.sphere_spectral_space = 'enable';
jg.compile.sphere_spectral_dealiasing = 'enable';

# Verbosity mode
jg.runtime.verbosity = 3

jg.runtime.output_file_mode = 'bin';

jg.runtime.benchmark_name = "galewsky"

#
# Compute error or difference to initial data
#
####jg.runtime.compute_error = 0

# Enable/Disbale GUI
jg = DisableGUI(jg)


jg.runtime.viscosity = 0

jg.runtime.reuse_plans = 'require_load'

max_simulation_time = 60 * 60 * 24 * 6.;

#
# Time, Mode and Physical resolution
#
timestep_size_reference = 60.;
timestep_size_fine = 60.; #3600 #1 hour  #864000/10 #1 day

jg.runtime.max_simulation_time = max_simulation_time; #1 day #timestep_size_reference #864000 #10 days
jg.runtime.output_timestep_size = max_simulation_time;
datastorage = jg.runtime.max_simulation_time / jg.runtime.output_timestep_size
if datastorage > 200:
	print("Warning::Too much data will be stored, are you sure you wish to run this?")

#jg.runtime.output_filename = "-"
#jg.runtime.output_timestep_size = timestep_size_reference*(2.0**(-timelevels))/10.0

jg.runtime.timestep_size = timestep_size_reference
jg.runtime.timestepping_method = "l_irk_n_erk"
jg.runtime.timestepping_order = 2
jg.runtime.timestepping_order2 = 2
jg.runtime.space_res_physical = -1
jg.runtime.space_res_spectral = 256

## Ref simulation
jg.reference_job = True
jg.compile.program = 'programs/pde_sweSphere'
jg.compile.parareal = "none";
jg.compile.xbraid = "none";
jg.runtime.parareal_enabled = 0;
jg.runtime.xbraid_enabled = 0;


params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
#params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
params_ptime_num_cores_per_rank = [1]

# Update TIME parallelization
ptime = JobParallelizationDimOptions('time')
ptime.num_cores_per_rank = 1
ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
ptime.num_ranks = 1

pspace = JobParallelizationDimOptions('space')
pspace.num_cores_per_rank = nb_cores_per_rank
###pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
pspace.num_threads_per_rank = nb_threads
pspace.num_ranks = 1

# Setup parallelization
jg.setup_parallelization([pspace, ptime])

jg.parallelization.mpiexec_disabled = False
####jg.parallelization.mpiexec_disabled = True


jg.gen_jobscript_directory();


## MGRIT simulations
jg.reference_job = False
jg.compile.program = "programs/xbraid_pde_sweSphere"
jg.reference_job_unique_id = jg.job_unique_id
ref_job = os.path.abspath(os.getcwd()) + "/" + jg.p_job_dirpath


jg.compile.xbraid = "mpi";
jg.runtime.xbraid_enabled = 1;
jg.runtime.xbraid_min_coarse = 2
jg.runtime.xbraid_nrelax0 = -1
jg.runtime.xbraid_tol = 0.
jg.runtime.xbraid_tnorm = 2
jg.runtime.xbraid_cfactor0 = -1
jg.runtime.xbraid_res = 0
jg.runtime.xbraid_storage = 0
jg.runtime.xbraid_print_level = 2
jg.runtime.xbraid_access_level = 2
jg.runtime.xbraid_run_wrapper_tests = 0
jg.runtime.xbraid_fullrnorm = 2
jg.runtime.xbraid_use_seq_soln = 0
jg.runtime.xbraid_use_rand = 1
jg.runtime.xbraid_timestepping_order = 2
jg.runtime.xbraid_timestepping_order2 = 2
jg.runtime.xbraid_verbosity = 0;

## path to ref job (to compute errors online!)
jg.runtime.xbraid_load_ref_csv_files = 1;
jg.runtime.xbraid_path_ref_csv_files = ref_job;
jg.runtime.xbraid_load_fine_csv_files = 1;
jg.runtime.xbraid_path_fine_csv_files = ref_job;
jg.runtime.xbraid_store_iterations = 1;

## Fixed relevant parameters
jg.runtime.timestep_size = timestep_size_reference;                           ## fine time step size
jg.runtime.xbraid_fmg = 1                                                     ## use F-cycles
jg.runtime.xbraid_fmg_vcyc = 1                                                ## use 1 post V-cycle after the F-cycle
jg.runtime.xbraid_skip = 1                                                    ## skip first down-cycle in order to compute the initial guess starting from the coarsest level
jg.runtime.xbraid_pt = 1                                                      ## number of parallel processors in time
jg.runtime.xbraid_max_iter = 11                                               ## number of iterations
jg.runtime.xbraid_no_output = 1                                               ## zero output!

cfactors = [2,4]; ## coarsening factor between levels
nbs_levels = [2,3] ## number of levels
nrelaxs = [0,1,2,5];
spatial_coarsening = [51,128];
tsms = ["l_irk_n_erk", "l_irk_na_sl_nr_settls_uv_only" ]


## test cases wall time
## tsm_coarse, nlevels, cfactor, nrelax,mcoarse,q1,nu1,q2,nu2
cases = [
            ["l_irk_na_sl_nr_settls_uv_only", 2, 2, 0, 51, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 2, 4, 0, 51, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 3, 2, 0, 51, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 3, 4, 0, 51, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 2, 2, 0, 128, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 2, 4, 0, 128, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 3, 2, 0, 128, 2, 1e7, 2, 1e7],
            ["l_irk_na_sl_nr_settls_uv_only", 3, 4, 0, 128, 2, 1e7, 2, 1e7],
        ]


for nb_pt in nb_pts:

    for [tsm_coarse, nb_levels, cfactor, nrelax, coarsening, q1, nu1, q2, nu2] in cases:

        jg.runtime.xbraid_timestepping_method = "l_irk_n_erk," + tsm_coarse;
        jg.runtime.xbraid_cfactor = cfactor;
        jg.runtime.xbraid_max_levels = nb_levels;
        jg.runtime.xbraid_pt = nb_pt;
        jg.runtime.xbraid_spatial_coarsening = coarsening;
        jg.runtime.xbraid_nrelax = nrelax;

        if nb_levels == 2:
            jg.runtime.xbraid_viscosity_order = "2," + str(q1);
            jg.runtime.xbraid_viscosity_coefficient = str(jg.runtime.viscosity) + "," + str(nu1);
        else:
            jg.runtime.xbraid_viscosity_order = "2," + str(q1) + "," + str(q2);
            jg.runtime.xbraid_viscosity_coefficient = str(jg.runtime.viscosity) + "," + str(nu1) + "," + str(nu2);

        params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
        #params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
        params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
        params_ptime_num_cores_per_rank = [1]

        # Update TIME parallelization
        ptime = JobParallelizationDimOptions('time')
        if nb_pt == 1:
            ptime.num_cores_per_rank = nb_cores_per_rank_pt1
            ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
            ptime.num_ranks = nb_pt
        else:
            ptime.num_cores_per_rank = 1
            ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
            ptime.num_ranks = nb_pt

        pspace = JobParallelizationDimOptions('space')
        pspace.num_cores_per_rank = nb_threads
        ###pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
        pspace.num_threads_per_rank = nb_threads
        pspace.num_ranks = 1

        # Setup parallelization
        jg.setup_parallelization([pspace, ptime])


        jg.gen_jobscript_directory()


