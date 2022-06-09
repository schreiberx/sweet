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
####tsm_fine = sys.argv[1];
####tsm_coarse = sys.argv[2];
simulation_to_run = sys.argv[3];
itest = int(sys.argv[4]);
###online_error = int(sys.argv[4])
###
###if online_error:
###    path_ref = sys.argv[5];
###    path_fine = sys.argv[6];

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

#
# Time, Mode and Physical resolution
#
jg.runtime.max_simulation_time = 1.
jg.runtime.output_timestep_size = .1

timestep_size_reference = 0.001
timestep_size_fine = 0.005
jg.runtime.timestep_size = timestep_size_fine
jg.runtime.space_res_spectral = 32


if simulation_to_run == "xbraid":

    jg.compile.xbraid = "mpi";
    jg.runtime.xbraid_enabled = 1;
    jg.runtime.xbraid_max_levels = 1
    jg.runtime.xbraid_skip = 1
    jg.runtime.xbraid_min_coarse = 1
    jg.runtime.xbraid_nrelax = 1
    jg.runtime.xbraid_nrelax0 = -1
    jg.runtime.xbraid_tol = 1e-9
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
    jg.runtime.xbraid_pt = 1
    jg.runtime.xbraid_timestepping_method = "dummy"
    jg.runtime.xbraid_timestepping_order = "2"
    jg.runtime.xbraid_timestepping_order2 = "2"


    if (itest == 0):
        jg.runtime.xbraid_run_wrapper_tests = 1;
    elif (itest == 1):
        jg.runtime.xbraid_max_levels = 1
    elif (itest == 2):
        jg.runtime.xbraid_max_levels = 1
        jg.runtime.xbraid_pt = 2
    elif (itest == 3):
        jg.runtime.xbraid_max_levels = 2
        jg.runtime.xbraid_tol = 0.
        jg.runtime.xbraid_max_iter = 3
        jg.runtime.xbraid_use_seq_soln = 1
    elif (itest == 4):
        jg.runtime.xbraid_print_level = 3
    elif (itest == 5):
        sys.exit("Not implemented!");



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
