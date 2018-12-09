#! /usr/bin/env python3

import os
import sys
import stat
import math

from SWEET import *
p = SWEETJobGeneration()


#
# Run simulation on plane or sphere
#
p.compile.program = 'swe_plane'

p.compile.plane_or_sphere = 'plane'

p.compile.plane_spectral_space = 'enable'
p.compile.plane_spectral_dealiasing = 'enable'
p.compile.sphere_spectral_space = 'disable'
p.compile.sphere_spectral_dealiasing = 'disable'

p.compile.numa_block_allocator = 0

# Verbosity mode
p.runtime.verbosity = 3

#
# Mode and Physical resolution
#
p.runtime.space_res_spectral = -1
p.runtime.space_res_physical = 512

#
# Benchmark ID
# 14: Steady diagonal benchmark
#
p.runtime.bench_id = 1

#
# Compute error
#
p.runtime.compute_error = 1

#
# Preallocate the REXI matrices
#
p.runtime.rexi_sphere_preallocation = 0

#
# Threading accross all REXI terms
#
rexi_thread_par = False
if rexi_thread_par:
	# OMP parallel for over REXI terms
	p.compile.threading = 'off'
	p.compile.rexi_thread_parallel_sum = 'enable'
else:
	p.compile.threading = 'omp'
	p.compile.rexi_thread_parallel_sum = 'disable'


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
# ==> use direct
p.runtime.rexi_method = 'terry'
#p.runtime.rexi_ci_n = 64
#p.runtime.rexi_ci_sx = 50
#p.runtime.rexi_ci_sy = 50
#p.runtime.rexi_ci_mu = 0
#p.runtime.rexi_ci_primitive = 'circle'
p.runtime.rexi_use_direct_solution = 1
         
# Parameters for SL-REXI paper
#-----------------------------       
p.runtime.gravitation= 9.80616
p.runtime.sphere_rotating_coriolis_omega = 0.00014584
p.runtime.h0 = 10000
p.runtime.plane_domain_size = 40031555.8928087

p.runtime.viscosity = 0.0

timelevels = 7 #5
timestep_size_reference = 3600 #1 hour  #864000/10 #1 day
timestep_sizes = [timestep_size_reference*(2.0**(-i)) for i in range(0, timelevels)]

p.runtime.max_simulation_time = 86400 #1 day #timestep_size_reference #864000 #10 days
p.runtime.output_timestep_size = p.runtime.max_simulation_time
#p.runtime.output_timestep_size = timestep_size_reference*(2.0**(-timelevels))/10.0

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
	prefix_string_template = group

	#
	# Reference solution
	#if True:
	if True:
		print("Reference")
		tsm = ts_methods[0]
	
		p.runtime.timestep_size = p.runtime.output_timestep_size/100.0
		p.runtime.timestepping_method = tsm[0]
		p.runtime.timestepping_order = tsm[1]
		p.runtime.timestepping_order2 = tsm[2]
		p.runtime.space_res_physical = 512

		if len(tsm) > 4:
			s = tsm[4]
			p.runtime.load_from_dict(tsm[4])

		p.gen_script('script_'+prefix_string_template+'_ref'+p.runtime.getUniqueID(p.compile), 'run.sh')

	for tsm in ts_methods[1:]:

		if group == 'ln2space' and 'ln_erk' in tsm[0]:
			p.runtime.space_grid_use_c_staggering = 1
			p.runtime.space_use_spectral_basis_diffs = 0

			#p.compile.plane_spectral_space = 'disable'
			p.compile.plane_spectral_dealiasing = 'disable'
			p.compile.libfft = 'enable'
		else:
			p.runtime.space_grid_use_c_staggering = 0
			p.runtime.space_use_spectral_basis_diffs = 1

			p.compile.plane_spectral_space = 'enable'
			p.compile.plane_spectral_dealiasing = 'enable'

		for idx in range(0, phys_res_levels): #, phys_res in phys_res_list:

			p.prefix_string = prefix_string_template

			p.runtime.timestep_size = timestep_sizes[idx]
			if group == 'ln2space' and 'ln_erk' in tsm[0]:
				p.runtime.timestep_size = p.runtime.timestep_size / 100.0

			p.runtime.timestepping_method = tsm[0]
			p.runtime.timestepping_order = tsm[1]
			p.runtime.timestepping_order2 = tsm[2]
			p.runtime.space_res_physical = phys_res_list[idx]
			print("id   dt       N  ")
			print(idx, p.runtime.timestep_size, p.runtime.space_res_physical)

			if len(tsm) > 4:
				s = tsm[4]
				p.runtime.load_from_dict(tsm[4])

			p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile), 'run.sh')

