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
from SWEETJobGeneration import *
from SWEETParameters import *

#Create main compile/run options
p = SWEETJobGeneration()

#Get Earth parameters (if necessary)
earth = EarthMKSDimensions()

#
# Run simulation on plane or sphere
#
#Basic plane options
p = CompileSWEPlane(p)


# Verbosity mode
p.runtime.verbosity = 3

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
#p.runtime.bench_id = 1
p.runtime.benchmark_name = "unstablejet"

#
# Compute error or difference to initial data
#
p.runtime.compute_error = 1

# Enable/Disbale GUI
#p = EnableGUI(p)
p = DisableGUI(p)

#
# REXI method
p.runtime.rexi_method = 'direct'
#p.runtime.rexi_use_direct_solution = 1

# Parameters for SL-REXI paper
#-----------------------------
p = RuntimeSWEPlaneEarthParam(p)
#p = RuntimeSWENonDimParam(p)


#
# Time, Mode and Physical resolution
#
timelevels = 6 #7 #5
timestep_size_reference = earth.day/24 #3600 #1 hour  #864000/10 #1 day
timestep_sizes = [timestep_size_reference*(2.0**(-i)) for i in range(0, timelevels)]

p.runtime.simtime = earth.day*12 #1 day #timestep_size_reference #864000 #10 days
p.runtime.output_timestep_size = p.runtime.simtime/24
datastorage = p.runtime.simtime / p.runtime.output_timestep_size
if datastorage > 200:
	print("Warning::Too much data will be stored, are you sure you wish to run this?")

#p.runtime.output_filename = "-"
#p.runtime.output_timestep_size = timestep_size_reference*(2.0**(-timelevels))/10.0

phys_res_levels = timelevels
phys_res_reference = 512
#phys_res_list = [phys_res_reference*(2**i) for i in range(0, phys_res_levels)]
phys_res_list = [phys_res_reference for i in range(0, phys_res_levels)]

p.runtime.viscosity = 0.0
p.runtime.viscosity_order = 2 #hyperviscosity
p.runtime.uselocalvisc = 1
visclevels = 10 #7 #5
visc_reference = 10000 #3600 #1 hour  #864000/10 #1 day
visc_sizes = [visc_reference*(10**(i)) for i in range(0, visclevels)]


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
			#['ln_erk',		2,	2],	# FD- C-grid
			#['l_cn_na_sl_nd_settls', 2,	2],	# SI-SL-SP
			['l_rexi_na_sl_nd_settls',	2,	2], #SL-EXP-SETTLS
			['l_rexi_na_sl_nd_etdrk',	2,	2], #SL-EXP-ETDRK
			#['l_rexi_n_etdrk',	2,	2], #ETDRK2
			#['l_rexi_n_erk',	2,	2], #strang split
		]

	#
	# OVERRIDE TS methods
	#
	if len(sys.argv) > 4:
		ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4])]]

	#
	# add prefix string to group benchmarks
	#
	prefix_string_template = group
	p.prefix_string = prefix_string_template

	#
	# Reference solution
	#if True:
	if False:
		print("Reference")
		tsm = ts_methods[0]
		
		p = SetupSpectralMethods(p)
		p.runtime.timestep_size = 2 # second #p.runtime.output_timestep_size/100.0
		p.runtime.timestepping_method = tsm[0]
		p.runtime.timestepping_order = tsm[1]
		p.runtime.timestepping_order2 = tsm[2]
		p.runtime.phys_res = -1
		p.runtime.mode_res = 1024

		if len(tsm) > 4:
			s = tsm[4]
			p.runtime.load_from_dict(tsm[4])

		p.gen_script('script_'+prefix_string_template+'_ref'+p.runtime.getUniqueID(p.compile), 'run.sh')

	for tsm in ts_methods[1:]:

		if group == 'sl-rexi' and 'ln_erk' in tsm[0]:
			p = SetupFDCMethods(p)
		else:
			p = SetupSpectralMethods(p)

		for idx in range(0, timelevels): #, phys_res in phys_res_list:

			p.runtime.timestep_size = timestep_sizes[idx]
			if group == 'sl-rexi' and 'ln_erk' in tsm[0]:
				p.runtime.timestep_size = p.runtime.timestep_size

			p.runtime.timestepping_method = tsm[0]
			p.runtime.timestepping_order = tsm[1]
			p.runtime.timestepping_order2 = tsm[2]
			p.runtime.phys_res = -1
			p.runtime.mode_res = phys_res_list[idx]
			print("id   dt       Nmodes  ")
			print(idx, p.runtime.timestep_size, p.runtime.phys_res)

			if len(tsm) > 4:
				s = tsm[4]
				p.runtime.load_from_dict(tsm[4])
				
			for ivisc in range(0, visclevels): 
				p.runtime.viscosity = visc_sizes[ivisc]
				p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile), 'run.sh')
