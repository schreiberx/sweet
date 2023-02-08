#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule.JobGeneration import *
jg = JobGeneration()


#
# Run simulation on plane or sphere
#
jg.compile.program = 'swe_plane'

jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'
jg.compile.sphere_spectral_space = 'disable'
jg.compile.sphere_spectral_dealiasing = 'disable'

jg.compile.numa_block_allocator = 0

# Verbosity mode
jg.runtime.verbosity = 3

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = -1
jg.runtime.space_res_physical = 512

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
jg.runtime.benchmark_name = "benchmark_id_1"
#jg.runtime.benchmark_name = "galewsky"

#
# Compute error
#
jg.runtime.compute_error = 1

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 0

#
# Threading accross all REXI terms
#
rexi_thread_par = False
if rexi_thread_par:
	# OMP parallel for over REXI terms
	jg.compile.threading = 'off'
	jg.compile.rexi_thread_parallel_sum = 'enable'
else:
	jg.compile.threading = 'omp'
	jg.compile.rexi_thread_parallel_sum = 'disable'


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
# ==> use direct
jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_ci_n = 64
#jg.runtime.rexi_ci_sx = 50
#jg.runtime.rexi_ci_sy = 50
#jg.runtime.rexi_ci_mu = 0
#jg.runtime.rexi_ci_primitive = 'circle'
#jg.runtime.rexi_use_direct_solution = 1
         
# Parameters for SL-REXI paper
#-----------------------------       
jg.runtime.gravitation= 9.80616
jg.runtime.sphere_rotating_coriolis_omega = 0.00014584
jg.runtime.h0 = 10000
jg.runtime.plane_domain_size = 40031555.8928087

jg.runtime.viscosity = 0.0

timelevels = 7 #5
timestep_size_reference = 3600 #1 hour  #864000/10 #1 day
timestep_sizes = [timestep_size_reference*(2.0**(-i)) for i in range(0, timelevels)]

jg.runtime.max_simulation_time = 86400 #1 day #timestep_size_reference #864000 #10 days
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
#jg.runtime.output_timestep_size = timestep_size_reference*(2.0**(-timelevels))/10.0

phys_res_levels = timelevels
phys_res_reference = 8
phys_res_list = [phys_res_reference*(2**i) for i in range(0, phys_res_levels)]


# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2']
groups = ['ln2space']


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:

	# 2nd order nonlinear non-fully-spectral
	if group == 'ln2space':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution - spectral (128 grid points)
			['ln_erk',		2,	2],	# FD- C-grid

			['l_cn_na_sl_nd_settls', 2,	2],	# SI-SL-SP

                        ['l_rexi_na_sl_nd_settls',	2,	2], #SL-EXP-SETTLS
			['l_rexi_na_sl_nd_etdrk',	2,	2], #SL-EXP-ETDRK

	#		['l_rexi_n_erk',	2,	2],
		]

	#
	# OVERRIDE TS methods
	#
	if len(sys.argv) > 4:
		ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4])]]

	#
	# add prefix string to group benchmarks
	#
	#prefix_string_template = group

	#
	# Reference solution
	#if True:
	if True:
		print("Reference")
		tsm = ts_methods[0]
	
		jg.runtime.timestep_size = min(timestep_sizes)
		jg.runtime.timestepping_method = tsm[0]
		jg.runtime.timestepping_order = tsm[1]
		jg.runtime.timestepping_order2 = tsm[2]
		jg.runtime.space_res_physical = 512

		if len(tsm) > 4:
			s = tsm[4]
			jg.runtime.load_from_dict(tsm[4])

		jg.reference_job = True
		jg.gen_jobscript_directory()
		jg.reference_job = False

		# Use this one as the reference solution!
		jg.reference_job_unique_id = jg.job_unique_id


	for tsm in ts_methods[1:]:

		if group == 'ln2space' and 'ln_erk' in tsm[0]:
			jg.runtime.space_grid_use_c_staggering = 1
			jg.runtime.space_use_spectral_basis_diffs = 0

			#jg.compile.plane_spectral_space = 'disable'
			jg.compile.plane_spectral_dealiasing = 'disable'
			jg.compile.libfft = 'enable'
		else:
			jg.runtime.space_grid_use_c_staggering = 0
			jg.runtime.space_use_spectral_basis_diffs = 1

			jg.compile.plane_spectral_space = 'enable'
			jg.compile.plane_spectral_dealiasing = 'enable'

		for idx in range(0, phys_res_levels): #, phys_res in phys_res_list:

			#jg.prefix_string = prefix_string_template

			jg.runtime.timestep_size = timestep_sizes[idx]
			if group == 'ln2space' and 'ln_erk' in tsm[0]:
				jg.runtime.timestep_size = jg.runtime.timestep_size / 100.0

			jg.runtime.timestepping_method = tsm[0]
			jg.runtime.timestepping_order = tsm[1]
			jg.runtime.timestepping_order2 = tsm[2]
			jg.runtime.space_res_physical = phys_res_list[idx]
			print("id   dt       N  ")
			print(idx, jg.runtime.timestep_size, jg.runtime.space_res_physical)

			if len(tsm) > 4:
				s = tsm[4]
				jg.runtime.load_from_dict(tsm[4])

			#jg.gen_script('script_'+prefix_string_template+jg.runtime.getUniqueID(jg.compile), 'run.sh')
			jg.gen_jobscript_directory()

