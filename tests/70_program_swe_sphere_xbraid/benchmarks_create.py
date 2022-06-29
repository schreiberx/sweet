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


####tsm_ref = "ln_erk";

simulation_to_run = sys.argv[1];
itest = int(sys.argv[2]);
tsm_fine = sys.argv[3];
tsm_coarse = sys.argv[4];

if (itest == 5):
    online_error = int(sys.argv[5])

    if online_error:
        path_fine = sys.argv[6];
###    path_ref = sys.argv[5];

#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
jg.compile.program = "swe_sphere"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"
jg.compile.fortran_source = "disable"

# Verbosity mode
jg.runtime.verbosity = 3

jg.compile.sphere_spectral_space = "enable";
jg.compile.sphere_spectral_dealiasing = "enable";

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "rossby_haurwitz_wave"

#
# Compute error or difference to initial data
#
jg.runtime.compute_error = 0

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
###jg.runtime.max_simulation_time = 500.
##jg.runtime.max_simulation_time = 3600.
jg.runtime.max_simulation_time = 720.
jg.runtime.output_timestep_size = 36.

timestep_size_reference = 18
timestep_size_fine = 36

jg.runtime.timestep_size = timestep_size_fine
jg.runtime.timestepping_method = tsm_fine
jg.runtime.timestepping_order = 2;
jg.runtime.timestepping_order2 = 2;


jg.runtime.space_res_spectral = 32

cfactors = [2, 4];
nbs_levels = [2, 4];
nb_pts = [1];
spatial_coarsening = [0, 1];

if simulation_to_run == "xbraid":

    jg.compile.xbraid = "mpi";
    jg.runtime.xbraid_enabled = 1;
    jg.runtime.xbraid_max_levels = 3
    jg.runtime.xbraid_skip = 1
    jg.runtime.xbraid_min_coarse = 2
    jg.runtime.xbraid_nrelax = 1
    jg.runtime.xbraid_nrelax0 = -1
    jg.runtime.xbraid_tol = 1e-9
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
    jg.runtime.xbraid_pt = 1
    jg.runtime.xbraid_timestepping_method = tsm_fine
    jg.runtime.xbraid_timestepping_order = 2
    jg.runtime.xbraid_timestepping_order2 = 2
    jg.runtime.xbraid_verbosity = 0;
    jg.runtime.xbraid_load_ref_csv_files = 0;
    jg.runtime.xbraid_path_ref_csv_files = "";
    jg.runtime.xbraid_load_fine_csv_files = 0;
    jg.runtime.xbraid_path_fine_csv_files = "";
    jg.runtime.xbraid_store_iterations = 0;
    jg.runtime.xbraid_spatial_coarsening = 0;

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
        jg.runtime.xbraid_timestepping_order = "2,2"
        jg.runtime.xbraid_timestepping_order2 = "2,2"
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
                    for coarsening in spatial_coarsening:
                        jg.runtime.xbraid_cfactor = cfactor;
                        jg.runtime.xbraid_max_levels = nb_levels;
                        jg.runtime.xbraid_pt = pt;
                        jg.runtime.xbraid_spatial_coarsening = coarsening;
                        jg.gen_jobscript_directory()

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
