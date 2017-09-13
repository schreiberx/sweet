#! /usr/bin/env python3

import os
import sys
import stat
import math

from SWEETJobGeneration import *
p = SWEETJobGeneration()


#
# Run simulation on plane or sphere
#
p.compile.program = 'swe_plane_rexi'

p.compile.plane_or_sphere = 'plane'

p.compile.plane_spectral_space = 'enable'
p.compile.plane_spectral_dealiasing = 'enable'
p.compile.sphere_spectral_space = 'disable'
p.compile.sphere_spectral_dealiasing = 'disable'

p.compile.numa_block_allocator = '2'

# Verbosity mode
p.runtime.verbosity = 3

#
# Mode and Physical resolution
#
p.runtime.mode_res = -1
p.runtime.phys_res = 512

#
# Benchmark ID
# 1: Gaussian breaking dam
#
p.runtime.bench_id = 14

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
rexi_thread_par = True
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
#
p.runtime.rexi_method = 'ci'
p.runtime.rexi_ci_n = 64
p.runtime.rexi_ci_sx = 50
p.runtime.rexi_ci_sy = 50
p.runtime.rexi_ci_mu = 0
p.runtime.rexi_ci_primitive = 'circle'


p.runtime.g = 1
p.runtime.f = 1
p.runtime.h = 100
p.runtime.domain_size = 1

p.runtime.viscosity = 0.0


timestep_size_reference = 0.0001
timestep_sizes = [0.0001*(2.0**i) for i in range(0, 11)]

p.runtime.simtime = 0.0001
p.runtime.output_timestep_size = p.runtime.simtime

phys_res_list = [16*(2**i) for i in range(0, 7)]

p.runtime.nonlinear = 1

# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2']
groups = ['ln2space']


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
	# 1st order linear
	if group == 'l1':
		ts_methods = [
			['l_direct',	0,	0,	0,	{'timestep_size': p.simtime}],	# reference solution
			['l_erk',	1,	0],
			['l_irk',	1,	0],
			['l_rexi',	0,	0],
		]

	# 2nd order linear
	if group == 'l2':
		ts_methods = [
			['l_direct',	0,	0,	{'timestep_size': p.simtime}],	# reference solution
			['l_erk',	2,	0],
			['l_cn',	2,	0],
			['l_rexi',	0,	0],
		]

	# 1st order nonlinear
	if group == 'ln1':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['l_erk_n_erk',		1,	1],
			['l_irk_n_erk',		1,	1],
			['ln_erk',		1,	1],
			['l_rexi_n_erk',	1,	1],
		]

	# 1st order nonlinear
	if group == 'ln1test':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['l_erk_n_erk',		1,	1],
			['l_irk_n_erk',		1,	1],
			['ln_erk',		1,	1],
		]

	# 2nd order nonlinear
	if group == 'ln2':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['l_cn_n_erk',		2,	2],
			['l_erk_n_erk',		2,	2],
			['l_irk_n_erk',		2,	2],
			['ln_erk',		2,	2],
			['l_rexi_n_erk',	2,	2],
		]

	# 2nd order nonlinear non-fully-spectral
	if group == 'ln2space':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution - spectral (128 grid points)
			['ln_erk',		2,	2],	# FD- C-grid
			['l_cn_na_sl_nd_settls', 2,	2],	# SI-SL-SP
#			['l_erk_n_erk',		2,	2],
#			['ln_erk',		2,	2],
#			['l_rexi_n_erk',	2,	2],
		]

	#
	# OVERRIDE TS methods
	#
	if len(sys.argv) > 4:
		ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]]


	#
	# add prefix string to group benchmarks
	#
	prefix_string_template = group


	#
	# Reference solution
	# (Not required, since we compare it with the real errors with --compute-errors=1)
	#if True:
	if False:
		print("Reference")
		tsm = ts_methods[0]
	
		p.timestep_size = timestep_size_reference
		p.runtime.timestepping_method = tsm[0]
		p.runtime.timestepping_order = tsm[1]
		p.runtime.timestepping_order2 = tsm[2]
		p.phys_res = 512

		if len(tsm) > 4:
			s = tsm[4]
			p.runtime.load_from_dict(tsm[4])

		p.gen_script('script_'+prefix_string_template+'_ref_'+p.runtime.getUniqueID(p.compile), 'run.sh')


	for tsm in ts_methods[1:]:

		if group == 'ln2space' and 'ln_erk' in tsm[0]:
			p.runtime.staggering = 1
			p.runtime.spectralderiv = 0

			p.compile.plane_spectral_space = 'disable'
			p.compile.plane_spectral_dealiasing = 'disable'
			p.compile.sphere_spectral_space = 'disable'
			p.compile.sphere_spectral_dealiasing = 'disable'
			p.compile.libfft = 'enable'

		if group == 'ln2space' and 'l_cn_na_sl_nd_settls' in tsm[0]:
			p.runtime.staggering = 0
			p.runtime.spectralderiv = 1

			p.compile.plane_spectral_space = 'enable'
			p.compile.plane_spectral_dealiasing = 'enable'
			p.compile.sphere_spectral_space = 'disable'
			p.compile.sphere_spectral_dealiasing = 'disable'

		for phys_res in phys_res_list:

			p.prefix_string = prefix_string_template

			p.runtime.timestep_size = timestep_size_reference
			p.runtime.timestepping_method = tsm[0]
			p.runtime.timestepping_order = tsm[1]
			p.runtime.timestepping_order2 = tsm[2]
			p.runtime.phys_res = phys_res


			if len(tsm) > 4:
				s = tsm[4]
				p.runtime.load_from_dict(tsm[4])

			p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile), 'run.sh')

