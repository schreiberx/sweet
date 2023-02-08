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


####tsm_ref = "ln_erk";
####tsm_fine = sys.argv[1];
####tsm_coarse = sys.argv[2];
simulation_to_run = sys.argv[3];
online_error = int(sys.argv[4])

if online_error:
    path_ref = sys.argv[5];
    path_fine = sys.argv[6];

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
###jg.compile.sweet_mpi = "enable"

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
###jg.runtime.timestepping_method = tsm_fine
###jg.runtime.timestepping_order = orders[tsm_fine]
###jg.runtime.timestepping_order2 = orders[tsm_fine]
###jg.runtime.space_res_physical = -1
jg.runtime.space_res_spectral = 32

## Parareal parameters
###jg.compile.parareal = "mpi";
jg.compile.parareal = "serial";
jg.runtime.parareal_enabled = 1
jg.runtime.parareal_convergence_threshold = -1
jg.runtime.parareal_verbosity = 6
jg.runtime.parareal_max_simulation_time = jg.runtime.max_simulation_time;
###jg.runtime.parareal_coarse_timestepping_order = orders[tsm_coarse];
###jg.runtime.parareal_coarse_timestepping_order2 = orders[tsm_coarse]

parareal_coarse_slices = [5, 10];
parareal_coarse_timesteps = [0.01, 0.05, -1]
##parareal_coarse_slices = [4, 6];
##parareal_coarse_timesteps = [15.,  30., -1]
##jg.runtime.parareal_coarse_timestepping_method = tsm_coarse;

if simulation_to_run == "parareal":
    if online_error:
        jg.runtime.parareal_store_iterations = 0;
        jg.runtime.parareal_load_ref_csv_files = 1;
        jg.runtime.parareal_path_ref_csv_files = path_ref;
        jg.runtime.parareal_load_fine_csv_files = 1;
        jg.runtime.parareal_path_fine_csv_files = path_fine;

    ## parareal simulations
    for nb_coarse_slices in parareal_coarse_slices:

        for coarse_timestep in parareal_coarse_timesteps:

            jg.runtime.parareal_coarse_slices = nb_coarse_slices;
            jg.runtime.parareal_coarse_timestep_size = coarse_timestep;
            jg.gen_jobscript_directory()

elif simulation_to_run == "ref":
    if not online_error:
        ## fine simulation
        jg.compile.parareal = "none";
        jg.runtime.parareal_enabled = 0;
        jg.gen_jobscript_directory();
        f = open("tmp_fine_sim.txt", "w");
        f.write(jg.job_dirpath);
        f.close();

        ## if ref simulation does not exist
        print (glob("job_benchref*"))
        if len(glob("job_benchref*")) == 0:
            ## ref simulation
            jg.compile.parareal = "none";
            jg.runtime.parareal_enabled = 0;
            ###jg.runtime.timestepping_method = tsm_ref
            jg.runtime.timestep_size = timestep_size_reference
            jg.reference_job = True;
            jg.gen_jobscript_directory();
            f = open("tmp_ref_sim.txt", "w");
            f.write(jg.job_dirpath);
            f.close();
        else:
            f = open("tmp_ref_sim.txt", "w");
            f.write(glob("job_benchref*")[0]);
            f.close();
