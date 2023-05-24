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

jg.compile.program = "programs/pde_sweSphere"
jg.compile.mode = "debug"
jg.compile.sweet_mpi = "enable"


jg.compile.sphere_spectral_space = 'enable';
jg.compile.sphere_spectral_dealiasing = 'enable';

# Verbosity mode
jg.runtime.verbosity = 3

jg.runtime.output_file_mode = 'bin';

#
# Benchmark ID
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "three_gaussian_bumps_phi_pint"

# Enable/Disbale GUI
jg = DisableGUI(jg)

#
# REXI method
jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1

jg.runtime.viscosity = 0

max_simulation_time = 60 * 2;

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

####### fine simulation
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
pspace.num_cores_per_rank = 1
###pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
pspace.num_threads_per_rank = 1
pspace.num_ranks = 1

# Setup parallelization
jg.setup_parallelization([pspace, ptime])

jg.parallelization.mpiexec_disabled = False
####jg.parallelization.mpiexec_disabled = True

jg.runtime.reuse_plans = 'save'

for dt in [60]:
    jg.runtime.timestep_size = dt

    jg.gen_jobscript_directory();

