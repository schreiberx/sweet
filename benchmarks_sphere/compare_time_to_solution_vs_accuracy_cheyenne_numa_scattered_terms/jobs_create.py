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


p.compile.compiler = 'intel'


#
# Use Intel MPI Compilers
#
p.compile.compiler_c_exec = 'mpicc'
p.compile.compiler_cpp_exec = 'mpicxx'
p.compile.compiler_fortran_exec = 'mpif90'


#
# Activate Fortran source
#
p.compile.fortran_source = 'enable'


#
# MPI?
#
p.compile.sweet_mpi = 'enable'


# Verbosity mode
p.runtime.verbosity = 2

#
# Mode and Physical resolution
#
p.runtime.mode_res = 128
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

if True:
	p.compile.threading = 'off'
	p.compile.rexi_thread_parallel_sum = 'disable'

else:
	#
	# WARNING: rexi_thread_par does not work yet!!!
	# MPI Ranks are clashing onthe same node with OpenMP Threads!
	#rexi_thread_par = True
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
#
p.runtime.rexi_method = 'ci'
p.runtime.rexi_ci_n = 64
p.runtime.rexi_ci_sx = 50
p.runtime.rexi_ci_sy = 50
p.runtime.rexi_ci_mu = 0
p.runtime.rexi_ci_primitive = 'circle'

#p.compile.debug_symbols = False


#p.runtime.g = 1
#p.runtime.f = 1
#p.runtime.h = 1
#p.runtime.domain_size = 1

p.runtime.viscosity = 0.0



timestep_size_reference = 100
#timestep_sizes = [timestep_size_reference*(2.0**i) for i in range(0, 11)]
timestep_sizes = [timestep_size_reference*(2**i) for i in range(0, 7)]

#p.runtime.simtime = timestep_sizes[-1]*10 #timestep_size_reference*2000
p.runtime.simtime = 100*(2**5)*10
p.runtime.output_timestep_size = p.runtime.simtime
#p.runtime.output_timestep_size = -1


# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
groups = ['ln2']


#
# MPI ranks
#
mpi_ranks = [2**i for i in range(0, 12+1)]



#
# allow including this file
#
if __name__ == "__main__":

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
				['l_rexi_n_etdrk',	1,	1,	0],
			]

		# 2nd order nonlinear
		if group == 'ln2':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution
				#['l_irk_n_erk',		2,	2,	0],
				#['l_cn_n_erk',		2,	2,	0],
				['l_erk_n_erk',		2,	2,	0],
				['l_irk_n_erk',		2,	2,	0],
				['ln_erk',		2,	2,	0],
				['l_rexi_n_etdrk',	2,	2,	0],
				['lg_irk_lc_n_erk',		2,	2,	0],
				['lg_rexi_lf_n_etdrk',	2,	2,	0],
			]

		# 4th order nonlinear
		if group == 'ln4':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution
				['l_rexi_n_etdrk',		4,	4,	0],
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
		# Parallelization models
		#
		# Use 18 cores for each MPI task even if only 1 thread is used
		# This avoid any bandwidth-related issues
		#
		# Parallelization model (18 threads per rank)
		p.cluster.pm_space_cores_per_mpi_rank = 18
		p.cluster.pm_time_cores_per_mpi_rank = 1


		#
		# Reference solution
		#
		if True:
		#if False:
			tsm = ts_methods[0]
			p.runtime.timestep_size  = timestep_sizes[0]

			p.runtime.timestepping_method = tsm[0]
			p.runtime.timestepping_order = tsm[1]
			p.runtime.timestepping_order2 = tsm[2]
			p.runtime.rexi_use_direct_solution = tsm[3]

			p.cluster.par_time_cores = 1

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

					c = 1
					#if True:
					if False:
						range_cores_single_socket = [1, 2, 4, 8, 12, 16, 18]
						range_cores_node = range_cores_single_socket + [18+i for i in range_cores_single_socket]
					else:
						range_cores_node = [18,36]

					range_cores = range_cores_node + [36*i for i in range(2, p.cluster.total_max_nodes)]

					if p.runtime.rexi_ci_n not in range_cores:
						range_cores.append(p.runtime.rexi_ci_n)
					range_cores.sort()

					if True:
						#for N in [64, 128]:
						for N in [64]:
							#for r in [25, 50, 75]:
							# Everything starting and above 40 results in significant errors
							#for r in [20, 30]:
							for r in [20]:
								p.runtime.load_from_dict({'rexi_method': 'ci', 'ci_n':N, 'ci_sx':r, 'ci_sy':r, 'half_poles':0})

								for p.cluster.par_time_cores in range_cores:
									if p.cluster.par_time_cores >= p.runtime.rexi_ci_n:
										# Generate only scripts with max number of cores
										p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')
										break

#					for p.cluster.par_time_cores in range_cores:
#						p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')

				else:
					p.cluster.par_time_cores = 1

					p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')



