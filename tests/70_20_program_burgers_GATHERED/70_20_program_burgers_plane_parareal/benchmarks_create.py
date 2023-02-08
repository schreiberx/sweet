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


tsm_ref = "ln_erk";
tsm_fine = sys.argv[1];
tsm_coarse = sys.argv[2];
simulation_to_run = sys.argv[3];
online_error = int(sys.argv[4])

if online_error:
    path_ref = sys.argv[5];
    path_fine = sys.argv[6];

orders = {};
####orders["l_irk_n_erk"] = 1;
####orders["l_irk"] = 1;
####orders["l_erk_n_erk"] = 2;
####orders["l_erk"] = 2;
####orders["l_cn_n_erk"] = 2;
####orders["l_cn"] = 2;
####orders["l_rexi_n_etdrk"] = 2;
####orders["l_rexi_n_erk"] = 2;
####orders["l_rexi"] = 2;
####orders["ln_erk"] = 2;
####orders["l_direct"] = 2;
####orders["l_rexi_na_sl_nd_settls"] = 2;
####orders["l_rexi_na_sl_nd_etdrk"] = 2;
####orders["l_cn_na_sl_nd_settls"] = 2;



#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options

#
# Run simulation on plane or sphere
#
jg.compile.program = 'burgers'
jg.compile.mode = "debug"

jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'
jg.compile.sphere_spectral_space = 'disable'
jg.compile.sphere_spectral_dealiasing = 'disable'

# Enable quad math per default for CI REXI method
jg.compile.quadmath = 'enable'

# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = None
jg.runtime.space_res_physical = 32

#jg.runtime.benchmark_name = "gaussian_bumps_phi_vrt_div"
jg.runtime.benchmark_name = "70"

# Viscosity
jg.runtime.viscosity = 0.01

# Simulation time
jg.runtime.max_simulation_time = 1.

# Output data
jg.runtime.output_timestep_size = 0.1

# Compute error
jg.runtime.compute_error = 0

# Enable/Disbale GUI
#jg = EnableGUI(jg)
jg = DisableGUI(jg)


timestep_size_reference = 0.001
timestep_size_fine = 0.005

jg.runtime.timestep_size = timestep_size_fine
jg.runtime.timestepping_method = tsm_fine
##jg.runtime.timestepping_order = orders[tsm_fine]
##jg.runtime.timestepping_order2 = orders[tsm_fine]
jg.runtime.space_res_physical = -1
jg.runtime.space_res_spectral = 32

## Parareal parameters
jg.compile.parareal = "serial";
jg.runtime.parareal_enabled = 1
jg.runtime.parareal_convergence_threshold = -1
jg.runtime.parareal_verbosity = 6
jg.runtime.parareal_max_simulation_time = jg.runtime.max_simulation_time;
##jg.runtime.parareal_coarse_timestepping_order = orders[tsm_coarse];
##jg.runtime.parareal_coarse_timestepping_order2 = orders[tsm_coarse]
jg.runtime.parareal_coarse_timestepping_order = 2;
jg.runtime.parareal_coarse_timestepping_order2 = 2;

parareal_coarse_slices = [5, 10];
parareal_coarse_timesteps = [0.01, 0.05]
##parareal_coarse_slices = [4, 6];
##parareal_coarse_timesteps = [15.,  30., -1]
jg.runtime.parareal_coarse_timestepping_method = tsm_coarse;

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
        f = open("fine_sim", "w");
        f.write(jg.job_dirpath);
        f.close();

        ## if ref simulation does not exist
        print (glob("job_benchref*"))
        if len(glob("job_benchref*")) == 0:
            ## ref simulation
            jg.compile.parareal = "none";
            jg.runtime.parareal_enabled = 0;
            jg.runtime.timestepping_method = tsm_ref
            jg.runtime.timestep_size = timestep_size_reference
            jg.reference_job = True;
            jg.gen_jobscript_directory();
            f = open("ref_sim", "w");
            f.write(jg.job_dirpath);
            f.close();
        else:
            f = open("ref_sim", "w");
            f.write(glob("job_benchref*")[0]);
            f.close();
