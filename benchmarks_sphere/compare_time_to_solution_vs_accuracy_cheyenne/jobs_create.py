#! /usr/bin/env python3

import os
import sys
import stat
import math

from SWEETJobGeneration import *
p = SWEETJobGeneration()

p.cluster.setupTargetMachine("cheyenne")


#
# Run simulation on plane or sphere
#
p.compile.program = 'swe_sphere_rexi'

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'



# Verbosity mode
p.runtime.verbosity = 2

#
# Mode and Physical resolution
#
p.runtime.mode_res = 64
p.runtime.phys_res = -1

#
# Benchmark ID
# 4: Gaussian breaking dam
#
p.runtime.bench_id = 4

#
# Compute error
#
p.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
p.runtime.rexi_sphere_preallocation = 1

#
# Deactivate stability checks
#
p.stability_checks = 0

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


#p.runtime.g = 1
#p.runtime.f = 1
#p.runtime.h = 1
#p.runtime.domain_size = 1

p.runtime.viscosity = 0.0


timestep_size_reference = 2
timestep_sizes = [timestep_size_reference*(2.0**i) for i in range(0, 11)]

p.runtime.simtime = timestep_size_reference*2000
p.runtime.output_timestep_size = p.runtime.simtime


# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']


#
# MPI ranks
#
mpi_ranks = [2**i for i in range(0, 12+1)]


####################################################
# WE FOCUS ON 2ND ORDER ACCURATE METHODS HERE
####################################################
groups = ['ln2']


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
	# 1st order linear
	if group == 'l1':
		ts_methods = [
			['l_erk',		4,	4,	0],	# reference solution
			['l_erk',	1,	0,	0],
			['l_irk',	1,	0,	0],
			['l_rexi',	0,	0,	0],
		]

		#if True:
		if False:
			for h in [0.1, 0.2, 0.3, 0.4, 0.5]:
			#for h in [0.1]:
				for M in [2**i for i in range(4, 11)]:
				#for M in [2**i for i in range(4, 5)]:
					ts_methods.append(['l_rexi',	0,	0,	0, {'rexi_method': 'terry', 'h':h, 'm':M}])

		#if True:
		if False:
			for testabs in [2**i for i in range(0, 3)]:
				for max_error in [1e-6, 1e-8, 1e-10, 1e-12]:
					ts_methods.append(['l_rexi',	0,	0,	0, {'rexi_method': 'file', 'file_test_abs':testabs, 'file_max_error':max_error}])


	# 2nd order linear
	if group == 'l2':
		ts_methods = [
			['l_erk',	4,	4,	0],	# reference solution
			['l_erk',	2,	0,	0],
			['l_cn',	2,	0,	0],
			['l_rexi',	0,	0,	0],
		]

	# 1st order nonlinear
	if group == 'ln1':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			['l_erk_n_erk',		1,	1,	0],
			['l_irk_n_erk',		1,	1,	0],
			['ln_erk',		1,	1,	0],
			['l_rexi_n_erk',	1,	1,	0],
			['ln_etdrk',		1,	1,	0],
		]

	# 2nd order nonlinear
	if group == 'ln2':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			#['l_irk_n_erk',		2,	2,	0],
			#['l_cn_n_erk',		2,	2,	0],
			#['l_erk_n_erk',		2,	2,	0],
			['ln_erk',		2,	2,	0],
#			['l_rexi_n_erk',	2,	2,	0],
			['ln_etdrk',		2,	2,	0],
		]

		if True:
			for N in [32, 64, 96, 128]:
				for r in [25, 50, 75, 100]:
					ts_methods.append(['ln_etdrk',	2,	2,	0, {'rexi_method': 'ci', 'ci_n':N, 'ci_sx':r, 'ci_sy':r}])

	# 4th order nonlinear
	if group == 'ln4':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			['ln_etdrk',		4,	4,	0],
			['ln_erk',		4,	4,	0],
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
	#
	if True:
		tsm = ts_methods[0]

		p.runtime.timestepping_method = tsm[0]
		p.runtime.timestepping_order = tsm[1]
		p.runtime.timestepping_order2 = tsm[2]
		p.runtime.rexi_use_direct_solution = tsm[3]

		if len(tsm) > 4:
			s = tsm[4]
			p.load_from_dict(tsm[4])

		p.gen_script('script_'+prefix_string_template+'_ref'+p.runtime.getUniqueID(p.compile), 'run.sh')


	#
	# Create job scripts
	#
	for tsm in ts_methods[1:]:
		for p.runtime.timestep_size in timestep_sizes:
#			for p.cluster.par_time_cores in mpi_ranks:
			p.runtime.timestepping_method = tsm[0]
			p.runtime.timestepping_order = tsm[1]
			p.runtime.timestepping_order2 = tsm[2]
			p.runtime.rexi_use_direct_solution = tsm[3]

			if len(tsm) > 4:
				s = tsm[4]
				p.runtime.load_from_dict(tsm[4])

			if 'etdrk' in p.runtime.timestepping_method:
				p.cluster.par_time_cores = p.runtime.rexi_ci_n
			else:
				p.cluster.par_time_cores = 1

			p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')



