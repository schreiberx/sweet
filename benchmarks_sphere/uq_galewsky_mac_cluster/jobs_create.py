#! /usr/bin/env python3

import os
import sys
import stat
import math
import itertools
import numpy as np

from SWEETJobGeneration import *
p = SWEETJobGeneration()

p.cluster.setupTargetMachine("mac-login-amd")

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


# Generate standard sparse grid points

def genDist(
	distribution_info,
	num_samples	# number of all samples
):
	if 'avg' in distribution_info:
		avg = distribution_info['avg']
		if 'perc' in distribution_info:
			perc = distribution_info['perc']
			return [avg-avg*perc, avg, avg+avg*perc]
		else:
			raise "Not supported"
	else:
		raise "Not supported"



def genDistInt(
	distribution_info,
	num_samples	# number of all samples
):
	return [int(x) for x in genDist(distribution_info, num_samples)]



######################################################################
# Sampling space: Simulation parameters
######################################################################

#
# Coriolis effect
#
avg_coriolis_omega = 7.292e-5
samples_coriolis_omega = genDist({'avg': avg_coriolis_omega, 'perc':0.1}, 1)

#
# Gravitation
#
avg_gravitation = 9.80616
samples_gravitation = genDist({'avg': avg_gravitation, 'perc':0.1}, 1)

#
# Average height
#
avg_h0 = 10000
samples_h0 = genDist({'avg': avg_h0, 'perc':0.1}, 1)

#
# Earth radius
#
avg_earth_radius = 6.37122e6
samples_earth_radius = genDist({'avg': avg_earth_radius, 'perc': 0.1}, 1)

#
# Viscosity
#
avg_viscosity = 1e5
samples_viscosity = genDist({'avg': avg_viscosity, 'perc': 0.1}, 1)

#
# Resolution
#
avg_mode_res = 128
samples_mode_res = genDistInt({'avg': avg_mode_res, 'perc': 0.1}, 1)

#
# Timestep size
#
avg_timestep_size = 60.0
samples_timestep_size = genDistInt({'avg': avg_timestep_size, 'perc': 0.1}, 1)



###################################
# Sampling space: Galewsky benchmark
###################################

#
# Average velocity
#
avg_galewsky_umax = 80.0
samples_galewsky_umax = genDist({'avg': avg_galewsky_umax, 'perc': 0.1}, 1)

#
# Amplification of bump
#
avg_galewsky_hamp = 120.0
samples_galewsky_hamp = genDist({'avg': avg_galewsky_hamp, 'perc': 0.1}, 1)

#
# Latitudinal placement of bump
#
avg_galewsky_phi2 = 0.25*math.pi
samples_galewsky_phi2 = genDist({'avg': avg_galewsky_phi2, 'perc': 0.1}, 1)



################ read SG points #####################
all_sg_points = []
with open("uq_code/grid_points.txt", "r") as file:
	grid_points = file.readlines()	
	
	for grid_point in grid_points:
		coordinate = grid_point.split()
		
		grid_points_dir = []
		for c in coordinate:
			grid_points_dir.append(float(c))

		all_sg_points.append(grid_points_dir)

file.close()

all_sg_points = np.array(all_sg_points)

samples_coriolis_omega 	= all_sg_points.T[0]
samples_gravitation 	= all_sg_points.T[1]
samples_earth_radius 	= all_sg_points.T[2]
samples_viscosity 	= all_sg_points.T[3]
samples_mode_res 	= all_sg_points.T[4]
samples_timestep_size 	= all_sg_points.T[5]
samples_galewsky_umax 	= all_sg_points.T[6]
samples_galewsky_hamp 	= all_sg_points.T[7]
samples_galewsky_phi2 	= all_sg_points.T[8]

print('samples_coriolis_omega')
print(samples_coriolis_omega)

print('samples_gravitation')
print(samples_gravitation)

print('samples_earth_radius')
print(samples_earth_radius)

print('samples_viscosity')
print(samples_viscosity)

print('samples_mode_res')
print(samples_mode_res)

print('samples_timestep_size')
print(samples_timestep_size)

print('samples_galewsky_umax')
print(samples_galewsky_umax)

print('samples_galewsky_hamp')
print(samples_galewsky_hamp)

print('samples_galewsky_phi2')
print(samples_galewsky_phi2)

#exit(0)
###################################################
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
# Benchmark name
#
# galewsky_nosetparam: Galewsky without overriding simulation parameters (gravity, avg. height, etc.)
p.runtime.benchmark_name = "galewsky_nosetparam"

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


timestep_size_reference = 10

#timestep_sizes = timestep_sizes[1:]
#print(timestep_sizes)
#sys.exit(1)

# SET time step size for implicit TS
p.runtime.timestep_size = 360


#p.runtime.simtime = timestep_sizes[-1]*10 #timestep_size_reference*2000
p.runtime.simtime = 432000 #timestep_size_reference*(2**6)*10
#p.runtime.output_timestep_size = p.runtime.simtime
#p.runtime.output_filename = ""
p.runtime.output_timestep_size = 432000

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

	ts_methods = [
		['ln_erk',		4,	4,	0],	# reference solution

		#
		# Choose implicit TS method for linear parts
		# with explicit RK for non-linear parts
		#
		['lg_irk_lc_n_erk_ver1',	2,	2,	0],
	]


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

		p.gen_script('script_ref'+p.runtime.getUniqueID(p.compile), 'run.sh')

		ts_methods = ts_methods[1:]


	#
	# Create job scripts
	#
	for tsm in ts_methods:
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

		p.runtime.timestepping_method = tsm[0]
		p.runtime.timestepping_order = tsm[1]
		p.runtime.timestepping_order2 = tsm[2]
		p.runtime.rexi_use_direct_solution = tsm[3]

		if len(tsm) > 4:
			s = tsm[4]
			p.runtime.load_from_dict(tsm[4])

		print("Number of parameter samples: "+str(len(all_sg_points)))

		for t in all_sg_points:
			p.runtime.f	 	= t[0]
			p.runtime.g 		= t[1]
			p.runtime.r 		= t[2]
			p.runtime.viscosity 	= t[3]
			p.runtime.mode_res 	= int(t[4])
			p.runtime.timestep_size 	= t[5]

			p.runtime.benchmark_galewsky_umax = t[6]
			p.runtime.benchmark_galewsky_hamp = t[7]
			p.runtime.benchmark_galewsky_phi2 = t[8]

			if not '_rexi' in p.runtime.timestepping_method:

				p.runtime.rexi_method = ''
				p.cluster.par_time_cores = 1
				p.gen_script('script'+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')

			else:

				# NOT YET USED, That's for Exponential Integrators
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
												p.gen_script('script'+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')
												break


