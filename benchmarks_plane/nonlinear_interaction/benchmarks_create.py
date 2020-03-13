#! /usr/bin/env python3
#
#  Create series of job to be run with sweet
#
#  Pedro Peixoto <pedrosp@ime.usp.br>
#  modified from Martin Schreiber initial job_create
#
#-------------------------------------------------------

#
# Usage on CoolMUC:
# Change the Script
# job_benchref_COMP_plspec_pldeal_numa2_fft_gnu_thomp_release_RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_ln_erk_tso4_tsob4_C2.000e+00_S864000_REXIDIR_M1024_N-001_X40031555.89280872_rob1_PAR_r00001_cpr028_tpr028_DIMS_dummy028_X
# manually to this:
#
#SBATCH --clusters=serial
#SBATCH --time=96:00:00
#
# if a job might take over 2 days to be executed

import os
import sys
import stat
import math

from sweet.SWEETRuntimeParametersScenarios import *
from mule_local.JobGeneration import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *

#Classes containing sweet compile/run basic option
#from mule_local.JobGeneration import *
#from sweet.SWEETRuntimeParametersScenarios import *

#Create main compile/run options
jg = JobGeneration()

# Request dedicated compile script
jg.compilecommand_in_jobscript = True

# Wallclock time
max_wallclock_seconds = 2*24*60*60
ref_max_wallclock_seconds = 48*60*60
jg.parallelization.max_wallclock_seconds = ref_max_wallclock_seconds

# HPC stuff
pspace = JobParallelizationDimOptions('space')
pspace.num_cores_per_rank =jg.platform_resources.num_cores_per_node/2
pspace.num_threads_per_rank = 8
pspace.num_ranks = 1
jg.setup_parallelization([pspace])


#Basic plane options
CompileSWEPlane(jg)

# Activate benchmark timers
jg.compile.benchmark_timings='enable'

# Verbosity mode
jg.runtime.verbosity = 3

# Benchmark ID
jg.runtime.benchmark_name = "normalmodes"

# Compute error or difference to initial data
jg.runtime.compute_error = 1

# Enable/Disbale GUI
EnableGUI(jg)
#DisableGUI(jg)

# REXI method
jg.runtime.rexi_method = 'ln_erk'

#Earth like parameters
earth = EarthMKSDimensionsApprox()
RuntimeSWEPlaneEarthParamApprox(jg)

#Setup time info
timestep_size_reference = earth.day/12 #3600 #1 hour  #864000/10 #1 day
jg.runtime.timestep_size = 60 #seconds
jg.runtime.max_simulation_time = earth.day*100 #1 day #timestep_size_reference #864000 #10 days
jg.runtime.output_timestep_size = jg.runtime.timestep_size #jg.runtime.max_simulation_time/100

#Output info
#jg.runtime.output_filename = "-" #do this to avoid dumping data files
jg.runtime.output_filename = "-"

datastorage = 0
if jg.runtime.output_timestep_size > 0:
	datastorage = jg.runtime.max_simulation_time / jg.runtime.output_timestep_size

if datastorage > 200 and len(jg.runtime.output_filename)>3 :
	print("Warning::Too much data will be stored, are you sure you wish to run this?")


#UNIQUE ID - what to cut out of naming!
unique_id_filter = []

# Compile
unique_id_filter.append('compile')

# Runtime
#unique_id_filter.append('runtime.disc_space')
#unique_id_filter.append('runtime.rexi')
#unique_id_filter.append('runtime.simparams')
#unique_id_filter.append('runtime.benchmark')

# Parallelization
unique_id_filter.append('parallelization')

jg.unique_id_filter = unique_id_filter

#Setup spectral space
jg.runtime.space_grid_use_c_staggering = 0
jg.runtime.space_use_spectral_basis_diffs = 1
jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'

#Setup method
jg.runtime.timestepping_method = 'ln_erk'
jg.runtime.timestepping_order = 4
jg.runtime.timestepping_order2 = 4
jg.runtime.space_res_physical = -1

#Setup  space info
jg.runtime.space_res_spectral = 64


#Viscosity
jg.runtime.viscosity = 0.0

#Banchmark to be used naming: waves_N_k0_k1_d0_dwest_deast_k0_k1_d0_dwest_deast_k0_k1_d0_dwest_deast
#jg.runtime.benchmark_normal_modes_case ="waves_1_3_2_0_100_0"
#jg.runtime.benchmark_normal_modes_case ="waves_2_1_3_0_100_0_7_2_100_0_0"
#jg.runtime.benchmark_normal_modes_case ="waves_1_1_3_100_0_0"
jg.runtime.benchmark_normal_modes_case ="waves_1_1_2_0_1_0"

# Tag this as a reference job
jg.reference_job = True
jg.gen_jobscript_directory()

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id

# Write compile script
jg.write_compilecommands()

