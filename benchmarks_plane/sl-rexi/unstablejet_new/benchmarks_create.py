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

#Classes containing sweet compile/run basic option
from mule_local.JobGeneration import *
from sweet.SWEETRuntimeParametersScenarios import *

#Create main compile/run options
jg = JobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
jg = CompileSWEPlane(jg)


# Verbosity mode
jg.runtime.verbosity = 3

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
timelevels = 10 #7 #5
timestep_size_reference = earth.day/12 #3600 #1 hour  #864000/10 #1 day
timestep_sizes = [timestep_size_reference*(2.0**(-i)) for i in range(0, timelevels)]

jg.runtime.max_simulation_time = earth.day*12 #1 day #timestep_size_reference #864000 #10 days
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time/24
datastorage = jg.runtime.max_simulation_time / jg.runtime.output_timestep_size
if datastorage > 200:
	print("Warning::Too much data will be stored, are you sure you wish to run this?")

#jg.runtime.output_filename = "-"
#jg.runtime.output_timestep_size = timestep_size_reference*(2.0**(-timelevels))/10.0

phys_res_levels = timelevels
phys_res_reference = 512
#phys_res_list = [phys_res_reference*(2**i) for i in range(0, phys_res_levels)]
phys_res_list = [phys_res_reference for i in range(0, phys_res_levels)]

# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2']
groups = ['sl-rexi']


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:

	# 2nd order nonlinear non-fully-spectral
	if group == 'sl-rexi':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['ln_erk',		2,	2],	# FD- C-grid
			['l_cn_na_sl_nd_settls', 2,	2],	# SI-SL-SP
			['l_rexi_na_sl_nd_settls',	2,	2], #SL-EXP-SETTLS
			['l_rexi_na_sl_nd_etdrk',	2,	2], #SL-EXP-ETDRK
			['l_rexi_n_etdrk',	2,	2], #ETDRK2
			['l_rexi_n_erk',	2,	2], #strang split
		]

	#
	# OVERRIDE TS methods
	#
	if len(sys.argv) > 4:
		ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4])]]

	#
	# Reference solution
	#if True:
	if False:
		print("Reference")
		tsm = ts_methods[0]
		
		p = SetupSpectralMethods(p)
		jg.runtime.timestep_size = 2 # second #jg.runtime.output_timestep_size/100.0
		jg.runtime.timestepping_method = tsm[0]
		jg.runtime.timestepping_order = tsm[1]
		jg.runtime.timestepping_order2 = tsm[2]
		jg.runtime.space_res_physical = -1
		jg.runtime.space_res_spectral = 1024

		if len(tsm) > 4:
			s = tsm[4]
			jg.runtime.load_from_dict(tsm[4])

		# Tag this as a reference job
		jg.reference_job = True
		jg.gen_jobscript_directory()
		jg.reference_job = False

		# Use this one as the reference solution!
		jg.reference_job_unique_id = jg.job_unique_id


	for tsm in ts_methods[1:]:

		if group == 'sl-rexi' and 'ln_erk' in tsm[0]:
			jg = SetupFDCMethods(jg)
		else:
			jg = SetupSpectralMethods(jg)

		for idx in range(0, timelevels): #, phys_res in phys_res_list:

			jg.runtime.timestep_size = timestep_sizes[idx]
			if group == 'sl-rexi' and 'ln_erk' in tsm[0]:
				jg.runtime.timestep_size = jg.runtime.timestep_size

			jg.runtime.timestepping_method = tsm[0]
			jg.runtime.timestepping_order = tsm[1]
			jg.runtime.timestepping_order2 = tsm[2]
			jg.runtime.space_res_physical = -1
			jg.runtime.space_res_spectral = phys_res_list[idx]
			print("id   dt       Nmodes  ")
			print(idx, jg.runtime.timestep_size, jg.runtime.space_res_physical)

			if len(tsm) > 4:
				s = tsm[4]
				jg.runtime.load_from_dict(tsm[4])

			jg.gen_jobscript_directory()
