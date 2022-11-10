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
from mule_local.JobGeneration import *
from mule_local.SWEETRuntimeParametersScenarios import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *

orders = {};
orders["l_irk_n_erk"] = 2;
orders["ln_erk"] = 2;
orders["l_cn_n_settls"] = 2;
orders["l_exp_n_etdrk"] = 2;


####tsm_ref = "ln_erk";

simulation_to_run = sys.argv[1];
model = sys.argv[2];
itest = int(sys.argv[3]);
tsm_fine = sys.argv[4];
tsm_coarse = sys.argv[5];
nb_pt = int(sys.argv[6]);

if (itest == 5 or itest == 6):
    online_error = int(sys.argv[7])

    if online_error:
        path_fine = sys.argv[8];

###    path_ref = sys.argv[5];

#Create main compile/run options
jg = JobGeneration()

####jg.compile.scalar_type = "complex"

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
jg.compile.program = "parareal_ode"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.runtime.ode_model = model;

# Verbosity mode
jg.runtime.verbosity = 3

jg.compile.sphere_spectral_space = "enable";
jg.compile.sphere_spectral_dealiasing = "enable";

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#jg.runtime.bench_id = 1
##jg.runtime.benchmark_name = "unstablejet"
if model == "ode1":
    jg.compile.scalar_type = "double"
    jg.compile.N_ode = 1
    jg.runtime.function_param_y0_real = 0.123
    jg.runtime.function_param_y0_imag = 0.
    jg.runtime.function_param_L = 1.
    jg.runtime.function_param_N = 1.
elif model == "SWE_triad":
    jg.compile.scalar_type = "complex"
    jg.compile.N_ode = 3
    y01 = 0.75
    y02 = 1.0
    y03 = 0.
    jg.runtime.function_param_y0_real = "{y01:.32f},{y02:.32f},{y03:.32f}".format(y01 = y01.real, y02 = y02.real, y03 = y03.real)
    jg.runtime.function_param_y0_imag = "{y01:.32f},{y02:.32f},{y03:.32f}".format(y01 = y01.imag, y02 = y02.imag, y03 = y03.imag)
    jg.runtime.function_param_L = "-2.3169272638598764,2.285662617093631,-0.05017777332080157"
    jg.runtime.function_param_N = "0.23644915578412212,-0.2359816957698708,0.005172082005927391"
    jg.runtime.function_param_extra = "-2.3169272638598764,2.285662617093631,-0.05017777332080157,0.24874903565537898" ## omega1, omega2, omega3, S123
else:
    sys.exit("Invalid model: " + model);


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

#
# Time, Mode and Physical resolution
#
jg.runtime.max_simulation_time = 1.
jg.runtime.output_timestep_size = .1

timestep_size_reference = 0.001
timestep_size_fine = 0.005
jg.runtime.timestep_size = timestep_size_fine
jg.runtime.timestepping_method = tsm_fine
jg.runtime.timestepping_order = orders[tsm_fine]
jg.runtime.timestepping_order2 = orders[tsm_fine]
jg.runtime.space_res_spectral = 32
cfactors = [2, 4, 8];
nbs_levels = [2, 4];
nb_pts = [1];

if simulation_to_run == "xbraid":

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
    jg.runtime.xbraid_timestepping_method = tsm_fine
    jg.runtime.xbraid_timestepping_order = "2"
    jg.runtime.xbraid_timestepping_order2 = "2"
    jg.runtime.xbraid_verbosity = 0;
    jg.runtime.xbraid_load_ref_csv_files = 0;
    jg.runtime.xbraid_path_ref_csv_files = "";
    jg.runtime.xbraid_load_fine_csv_files = 0;
    jg.runtime.xbraid_path_fine_csv_files = "";
    jg.runtime.xbraid_store_iterations = 0;


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
        jg.setup_parallelization([pspace, ptime])


    if (itest == 0):
        jg.runtime.xbraid_run_wrapper_tests = 1;
    elif (itest == 1):
        jg.runtime.xbraid_max_levels = 1
        jg.runtime.xbraid_store_iterations = 1;
    elif (itest == 2):
        jg.runtime.xbraid_max_levels = 1
        jg.runtime.xbraid_store_iterations = 1;
        jg.runtime.xbraid_pt = 2
    elif (itest == 3):
        jg.runtime.xbraid_max_levels = 2
        jg.runtime.xbraid_tol = 0.
        jg.runtime.xbraid_max_iter = 3
        jg.runtime.xbraid_use_seq_soln = 1
        jg.runtime.xbraid_access_level = 2
    elif (itest == 4):
        jg.runtime.xbraid_print_level = 3
    elif (itest == 5):
        jg.runtime.xbraid_timestepping_method = tsm_fine + "," + tsm_coarse;
        jg.runtime.xbraid_store_iterations = 1;
        jg.runtime.xbraid_access_level = 2;

        if online_error:
            jg.runtime.xbraid_store_iterations = 0;
            ###jg.runtime.xbraid_load_ref_csv_files = 1;
            ###jg.runtime.xbraid_path_ref_csv_files = path_ref;
            jg.runtime.xbraid_load_fine_csv_files = 1;
            jg.runtime.xbraid_path_fine_csv_files = path_fine;

        ## parareal simulations
        for cfactor in cfactors:
            for nb_levels in nbs_levels:
                for pt in nb_pts:
                    jg.runtime.xbraid_cfactor = cfactor;
                    jg.runtime.xbraid_max_levels = nb_levels;
                    jg.runtime.xbraid_pt = pt;
                    jg.gen_jobscript_directory()

    ## xbraid vs parareal
    elif (itest == 6):

        for cfactor in cfactors:
            ## xbraid with two leves and F-relaxation only
            jg.compile.xbraid = "mpi";
            jg.runtime.xbraid_enabled = 1;
            jg.compile.parareal = "none";
            jg.runtime.parareal_enabled = 0;
            jg.runtime.parareal_max_iter = jg.runtime.xbraid_max_iter;
            jg.runtime.xbraid_access_level = 2;
            jg.runtime.xbraid_store_iterations = 1;
            jg.runtime.xbraid_load_fine_csv_files = 1;
            jg.runtime.xbraid_path_fine_csv_files = path_fine;
            jg.runtime.xbraid_max_levels = 2;
            jg.runtime.xbraid_skip = 1;
            jg.runtime.xbraid_fmg = 0;
            jg.runtime.xbraid_nrelax = 0; ## F-relaxation
            jg.runtime.xbraid_cfactor = cfactor;
            jg.gen_jobscript_directory();

            ## parareal
            jg.compile.xbraid = "none";
            jg.runtime.xbraid_enabled = 0;
            jg.compile.parareal = "mpi";
            jg.runtime.parareal_enabled = 1;
            jg.runtime.parareal_max_iter = jg.runtime.xbraid_max_iter;
            jg.runtime.parareal_coarse_timestepping_method = tsm_coarse
            jg.runtime.parareal_coarse_timestepping_order = 2
            jg.runtime.parareal_coarse_timestepping_order2 = 2
            jg.runtime.parareal_load_fine_csv_files = 1;
            jg.runtime.parareal_path_fine_csv_files = path_fine;
            jg.runtime.parareal_store_iterations = 1;
            jg.runtime.parareal_coarse_timestep_size = -1
            jg.runtime.parareal_coarse_slices = int(jg.runtime.max_simulation_time / (timestep_size_fine * cfactor) );
            jg.gen_jobscript_directory();

    if (itest < 5):
        jg.gen_jobscript_directory();

elif simulation_to_run == "ref":

    jg.compile.xbraid = "none";
    jg.runtime.xbraid_enabled = 0;

    ###if not online_error:
    ## fine simulation
    jg.compile.parareal = "none";
    jg.runtime.parareal_enabled = 0;
    jg.gen_jobscript_directory();
    f = open("fine_sim", "w");
    f.write(jg.job_dirpath);
    f.close();

    ####    ## if ref simulation does not exist
    ####    print (glob("job_benchref*"))
    ####    if len(glob("job_benchref*")) == 0:
    ####        ## ref simulation
    ####        jg.compile.parareal = "none";
    ####        jg.runtime.parareal_enabled = 0;
    ####        ###jg.runtime.timestepping_method = tsm_ref
    ####        jg.runtime.timestep_size = timestep_size_reference
    ####        jg.reference_job = True;
    ####        jg.gen_jobscript_directory();
    ####        f = open("ref_sim", "w");
    ####        f.write(jg.job_dirpath);
    ####        f.close();
    ####    else:
    ####        f = open("ref_sim", "w");
    ####        f.write(glob("job_benchref*")[0]);
    ####        f.close();
