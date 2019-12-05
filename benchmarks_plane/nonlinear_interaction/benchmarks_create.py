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

#Classes containing sweet compile/run basic option
from mule_local.JobGeneration import *
from sweet.SWEETRuntimeParametersScenarios import *

#Create main compile/run options
jg = JobGeneration()

# Request dedicated compile script
jg.compilecommand_in_jobscript = True

# Wallclock time
max_wallclock_seconds = 2*24*60*60
ref_max_wallclock_seconds = 48*60*60

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
CompileSWEPlane(jg)


# Activate benchmark timers
jg.compile.benchmark_timings='enable'



# Verbosity mode
jg.runtime.verbosity = 3

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#jg.runtime.bench_id = 1
jg.runtime.benchmark_name = "normalmodes"

#
# Compute error or difference to initial data
#
jg.runtime.compute_error = 1

# Enable/Disbale GUI
EnableGUI(jg)
#DisableGUI(jg)

#
# REXI method
jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_use_direct_solution = 1

# Parameters for SL-REXI paper
#-----------------------------
RuntimeSWEPlaneEarthParam(jg)
#RuntimeSWENonDimParam(jg)

jg.runtime.viscosity = 0.0


#
# Time, Mode and Physical resolution
#
timelevels = 1 #7 #5
timestep_size_reference = earth.day/12 #3600 #1 hour  #864000/10 #1 day
timestep_sizes = [timestep_size_reference*(2.0**(-i)) for i in range(0, timelevels)]

# Use 10 days as the reference solution
jg.runtime.max_simulation_time = earth.day*10 #1 day #timestep_size_reference #864000 #10 days
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
#datastorage = jg.runtime.max_simulation_time / jg.runtime.output_timestep_size
#if datastorage > 200:
#	print("Warning::Too much data will be stored, are you sure you wish to run this?")

#jg.runtime.output_filename = "-"
#jg.runtime.output_timestep_size = timestep_size_reference*(2.0**(-timelevels))/10.0

phys_res_levels = timelevels
phys_res_reference = 128 #512
#phys_res_list = [phys_res_reference*(2**i) for i in range(0, phys_res_levels)]
phys_res_list = [phys_res_reference for i in range(0, phys_res_levels)]

jg.runtime.benchmark_normal_modes_case ="single_2_1_1_0_0"

ts_methods = [
	['ln_erk',		4,	4]#,	# reference solution
	#['ln_erk',		2,	2],	# FD- C-grid
	#['l_cn_na_sl_nd_settls', 2,	2],	# SI-SL-SP
	#['l_rexi_na_sl_nd_settls',	2,	2], #SL-EXP-SETTLS
	#['l_rexi_na_sl_nd_etdrk',	2,	2], #SL-EXP-ETDRK
	#['l_rexi_n_etdrk',	2,	2], #ETDRK2
	#['l_rexi_n_erk',	2,	2], #strang split
]


unique_id_filter = []

# Compile
unique_id_filter.append('compile')

# Runtime
unique_id_filter.append('runtime.disc_space')
unique_id_filter.append('runtime.rexi')
unique_id_filter.append('runtime.simparams')
unique_id_filter.append('runtime.benchmark')

# Parallelization
unique_id_filter.append('parallelization')

jg.unique_id_filter = unique_id_filter

#
# Reference solution
if True:
#if False:
	print("Reference")
	tsm = ts_methods[0]
	
	jg.parallelization.max_wallclock_seconds = ref_max_wallclock_seconds

	SetupSpectralMethods(jg)
	jg.runtime.timestep_size = 2 # second #jg.runtime.output_timestep_size/100.0
	jg.runtime.timestepping_method = tsm[0]
	jg.runtime.timestepping_order = tsm[1]
	jg.runtime.timestepping_order2 = tsm[2]
	jg.runtime.space_res_physical = -1
	#jg.runtime.space_res_spectral = 1024
	jg.runtime.space_res_spectral = phys_res_list[0]

	# Tag this as a reference job
	jg.reference_job = True
	jg.gen_jobscript_directory()
	jg.reference_job = False

	# Use this one as the reference solution!
	jg.reference_job_unique_id = jg.job_unique_id



#
# Use only 2 iterations for Semi-Lagrangian methods
#
unique_id_filter.append('runtime.semi_lagrangian')
jg.runtime.semi_lagrangian_iterations = 2
jg.runtime.semi_lagrangian_convergence_threshold = -1


jg.parallelization.max_wallclock_seconds = max_wallclock_seconds

for tsm in ts_methods[1:]:

	if 'ln_erk' in tsm[0]:
		SetupFDCMethods(jg)
	else:
		SetupSpectralMethods(jg)

	for idx in range(0, timelevels): #, phys_res in phys_res_list:

		jg.runtime.timestep_size = timestep_sizes[idx]

		jg.runtime.timestepping_method = tsm[0]
		jg.runtime.timestepping_order = tsm[1]
		jg.runtime.timestepping_order2 = tsm[2]
		jg.runtime.space_res_physical = -1
		jg.runtime.space_res_spectral = phys_res_list[idx]
		print("id   dt       Nmodes  ")
		print(idx, jg.runtime.timestep_size, jg.runtime.space_res_physical)

		jg.gen_jobscript_directory()



# Write compile script
jg.write_compilecommands()

