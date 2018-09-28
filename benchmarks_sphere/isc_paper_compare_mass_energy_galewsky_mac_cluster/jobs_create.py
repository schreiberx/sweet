#! /usr/bin/env python3

import os
import sys
import stat
import math

sys.path.append(os.environ['SWEET_ROOT']+'/python_mods/')
from SWEETJobGeneration import *
p = SWEETJobGeneration()

p.cluster.setupTargetMachine("mac-login-amd")


#
# Run simulation on plane or sphere
#
p.compile.program = 'swe_sphere'

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'


#p.compile.compiler = 'intel'


#
# Use Intel MPI Compilers
#
#p.compile.compiler_c_exec = 'mpicc'
#p.compile.compiler_cpp_exec = 'mpicxx'
#p.compile.compiler_fortran_exec = 'mpif90'


#
# Activate Fortran source
#
p.compile.fortran_source = 'enable'


#
# MPI?
#
#p.compile.sweet_mpi = 'enable'


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
# 100: Galewski
#
p.runtime.bench_id = 100

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
	#p.compile.rexi_thread_parallel_sum = 'disable'
	p.compile.rexi_thread_parallel_sum = 'enable'

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
p.runtime.rexi_method = ''
p.runtime.rexi_ci_n = 128
p.runtime.rexi_ci_max_real = -999
p.runtime.rexi_ci_max_imag = -999
p.runtime.rexi_ci_sx = -1
p.runtime.rexi_ci_sy = -1
p.runtime.rexi_ci_mu = 0
p.runtime.rexi_ci_primitive = 'circle'

#p.runtime.rexi_beta_cutoff = 1e-16
p.runtime.rexi_beta_cutoff = 0

#p.compile.debug_symbols = False


#p.runtime.g = 1
#p.runtime.f = 1
#p.runtime.h = 1
#p.runtime.domain_size = 1

p.runtime.viscosity = 0.0



#timestep_sizes = [timestep_size_reference*(2.0**i) for i in range(0, 11)]
#timestep_sizes = [timestep_size_reference*(2**i) for i in range(2, 4)]

timestep_sizes_explicit = [10, 20, 30, 60, 120, 180]
timestep_sizes_implicit = [60, 120, 180, 360, 480, 600, 720]
timestep_sizes_rexi = [60, 120, 180, 240, 300, 360, 480, 600, 720]

timestep_size_reference = timestep_sizes_explicit[0]

#timestep_sizes = timestep_sizes[1:]
#print(timestep_sizes)
#sys.exit(1)


#p.runtime.simtime = timestep_sizes[-1]*10 #timestep_size_reference*2000
p.runtime.simtime = 432000 #timestep_size_reference*(2**6)*10
#p.runtime.output_timestep_size = p.runtime.simtime
p.runtime.output_filename = "-"
p.runtime.output_timestep_size = 60

p.runtime.rexi_extended_modes = 0

p.runtime.floating_point_output_digits = 14

# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
groups = ['ln2']


#
# MPI ranks
#
#mpi_ranks = [2**i for i in range(0, 12+1)]
#mpi_ranks = [1]



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

		# 2nd order nonlinear
		if group == 'ln2':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution

				###########
				# RK2/4
				###########
				['ln_erk',		2,	2,	0],	# reference solution
				['ln_erk',		4,	4,	0],	# reference solution

				###########
				# CN
				###########
				['lg_irk_lc_n_erk_ver0',	2,	2,	0],
				['lg_irk_lc_n_erk_ver1',	2,	2,	0],

				['l_irk_n_erk_ver0',	2,	2,	0],
				['l_irk_n_erk_ver1',	2,	2,	0],

				###########
				# REXI
				###########
				['lg_rexi_lc_n_erk_ver0',	2,	2,	0],
				['lg_rexi_lc_n_erk_ver1',	2,	2,	0],

				['l_rexi_n_erk_ver0',	2,	2,	0],
				['l_rexi_n_erk_ver1',	2,	2,	0],

				###########
				# ETDRK
				###########
				['lg_rexi_lc_n_etdrk',	2,	2,	0],
				['l_rexi_n_etdrk',	2,	2,	0],
			]

		# 4th order nonlinear
		if group == 'ln4':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution
				['l_rexi_n_etdrk',	4,	4,	0],
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
#		p.cluster.pm_space_cores_per_mpi_rank = 18
#		p.cluster.pm_time_cores_per_mpi_rank = 1


		#
		# Reference solution
		#
		if True:
		#if False:
			tsm = ts_methods[0]
			p.runtime.timestep_size  = timestep_size_reference

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
			tsm_name = tsm[0]
			if 'ln_erk' in tsm_name:
				timestep_sizes = timestep_sizes_explicit
			elif 'l_erk' in tsm_name or 'lg_erk' in tsm_name:
				timestep_sizes = timestep_sizes_explicit
			elif 'l_irk' in tsm_name or 'lg_irk' in tsm_name:
				timestep_sizes = timestep_sizes_implicit
			elif 'l_rexi' in tsm_name or 'lg_rexi' in tsm_name:
				timestep_sizes = timestep_sizes_rexi
			else:
				print("Unable to identify time stepping method "+tsm_name)
				sys.exit(1)

			for p.runtime.timestep_size in timestep_sizes:
				p.runtime.timestepping_method = tsm[0]
				p.runtime.timestepping_order = tsm[1]
				p.runtime.timestepping_order2 = tsm[2]
				p.runtime.rexi_use_direct_solution = tsm[3]

				if len(tsm) > 4:
					s = tsm[4]
					p.runtime.load_from_dict(tsm[4])

				if not '_rexi' in p.runtime.timestepping_method:
					p.runtime.rexi_method = ''
					p.cluster.par_time_cores = 1
					p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')

				else:
					c = 1
					#if True:
					if False:
						range_cores_single_socket = [1, 2, 4, 8, 12, 16, 18]
						range_cores_node = range_cores_single_socket + [18+i for i in range_cores_single_socket]
					else:
						#range_cores_node = [18,36]
						range_cores_node = [18]

					if True:
						#for N in [64, 128]:
						#for N in [128, 256]:
						#for N in [128, 256]:
						for N in [128]:

							range_cores = range_cores_node + [36*i for i in range(2, p.cluster.total_max_nodes)]

							if p.runtime.rexi_ci_n not in range_cores:
								range_cores.append(N)
							range_cores.sort()

							#for r in [25, 50, 75]:
							# Everything starting and above 40 results in significant errors
							#for r in [30, 50]:
							#for r in [30, 60]:
							for r in [30]:
								#for gf in [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.0]:
								#for gf in [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.0]:
								#for gf in [1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001]:
								#for gf_exp_N in [2, 4, 6, 10, 20, 40]:
								#for gf_exp_N in [2, 4, 10]:
								for gf_exp_N in [0]:
									#for gf_scale in [0, 5, 10, 20, 50]:
									for gf_scale in [0]:

										#for ci_max_real in [10, 5]:
										for ci_max_real in [10.0]:
											p.runtime.load_from_dict({
												'rexi_method': 'ci',
												'ci_n':N,
												'ci_max_real':ci_max_real,
												'ci_max_imag':r,
												'half_poles':0,
												'ci_gaussian_filter_scale':gf_scale,
												#'ci_gaussian_filter_dt_norm':130.0,	# unit scaling for T128 resolution
												'ci_gaussian_filter_dt_norm':0.0,	# unit scaling for T128 resolution
												'ci_gaussian_filter_exp_N':gf_exp_N,
											})

											#print(range_cores)
											#sys.exit(1)
											range_cores = [1]
											for p.cluster.par_time_cores in range_cores:
												#if p.cluster.par_time_cores >= p.runtime.rexi_ci_n:
												if True:
													# Generate only scripts with max number of cores
													p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')
													break

#					for p.cluster.par_time_cores in range_cores:
#						p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')


