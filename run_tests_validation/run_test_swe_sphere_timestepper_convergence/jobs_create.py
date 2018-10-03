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
p.compile.program = 'swe_sphere'

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'


# Enable quad math per default for CI REXI method
p.compile.quadmath = 'enable'


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
p.runtime.rexi_ci_n = 128
p.runtime.rexi_ci_max_real = 10
p.runtime.rexi_ci_max_imag = 20
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
groups = ['l1', 'lg1', 'l2', 'lg2', 'ln1', 'ln2', 'ln4']

if len(sys.argv) > 1:
	groups = [sys.argv[1]]


#
# allow including this file
#
if __name__ == "__main__":

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

			if False:
				ts_methods.append(['l_erk',	4,	0,	0,	{'timestep_size': 0.01}])

				for h in [0.1, 0.2, 0.3, 0.4, 0.5]:
				#for h in [0.1]:
					for M in [2**i for i in range(4, 11)]:
					#for M in [2**i for i in range(4, 5)]:
						ts_methods.append(['l_rexi',	0,	0,	0, {'rexi_method': 'terry', 'h':h, 'm':M}])

				for testabs in [2**i for i in range(0, 3)]:
					for max_error in [1e-6, 1e-8, 1e-10, 1e-12]:
						ts_methods.append(['l_rexi',	0,	0,	0, {'rexi_method': 'file', 'file_test_abs':testabs, 'file_max_error':max_error}])

		# 1st order linear
		if group == 'lg1':
			ts_methods = [
				['lg_erk',		4,	4,	0],	# reference solution
				['lg_erk',	1,	0,	0],
				['lg_irk',	1,	0,	0],
				#['l_rexi',	0,	0,	0],
			]

		# 2nd order linear
		if group == 'l2':
			ts_methods = [
				['l_erk',		4,	4,	0],	# reference solution
				['l_erk',	2,	0,	0],
				['l_cn',	2,	0,	0],
				['l_rexi',	0,	0,	0],
			]

		# 2nd order linear
		if group == 'lg2':
			ts_methods = [
				['lg_erk',		4,	4,	0],	# reference solution
				['lg_erk',	2,	0,	0],
				['lg_cn',	2,	0,	0],
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
				['l_rexi_n_etdrk',	1,	1,	0],
				['lg_rexi_lf_n_etdrk',	1,	1,	0],
			]

		# 2nd order nonlinear
		if group == 'ln2':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution
				['l_cn_n_erk',		2,	2,	0],
				['l_erk_n_erk',		2,	2,	0],
				['l_irk_n_erk',		2,	2,	0],
				['ln_erk',		2,	2,	0],
				['l_rexi_n_erk',	2,	2,	0],
				['l_rexi_n_etdrk',	2,	2,	0],
				['lg_rexi_lf_n_etdrk',	2,	2,	0],
			]

		# 4th order nonlinear
		if group == 'ln4':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution
				['l_rexi_n_etdrk',	4,	4,	0],
				['lg_rexi_lf_n_etdrk',	4,	4,	0],
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

			p.runtime.timestep_size = timestep_size_reference
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
				p.runtime.timestepping_method = tsm[0]
				p.runtime.timestepping_order = tsm[1]
				p.runtime.timestepping_order2 = tsm[2]
				p.runtime.rexi_use_direct_solution = tsm[3]

				if len(tsm) > 4:
					s = tsm[4]
					p.runtime.load_from_dict(tsm[4])

				p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile), 'run.sh')

