#! /usr/bin/env python3

import sys
import math

from mule import *
p = JobGeneration()


#
# Cluster options
#
p.cluster.setupTargetMachine("cheyenne")
p.cluster.pm_space_cores_per_mpi_rank = 18
p.cluster.pm_time_cores_per_mpi_rank = 1

par_time_cores_list = [2**i for i in range(0, 10)]



#
# Compile options
#


#
# PDE
#
p.compile.program = 'swe_sphere'

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'


#
# Use Intel MPI Compilers
#
p.compile.compiler_c_exec = 'mpicc'
p.compile.compiler_cpp_exec = 'mpicxx'
p.compile.compiler_fortran_exec = 'mpif90'


#
# parallelization
#
#p.compile.rexi_thread_parallel_sum = 'enable'
p.compile.rexi_thread_parallel_sum = 'disable'
p.compile.threading = 'off'
p.compile.sweet_mpi = 'enable'


p.compile.compiler = 'intel'
p.compile.fortran_source = 'enable'



p.runtime.timestep_size = 0.001
p.runtime.space_res_spectral = 128
p.runtime.output_timestep_size = 129600
p.runtime.max_simulation_time = 129600
#p.runtime.output_filename = '-'
p.runtime.timestepping_order = 4

p.runtime.timestepping_method = 'l_erk'
p.runtime.rexi_m = 256
p.runtime.rexi_h = 0.15
p.runtime.rexi_half_poles = 0
p.runtime.sphere_extended_modes = 0
p.runtime.rexi_normalization = 1

p.runtime.f_sphere = 0

p.runtime.gravitation= 9.80616	# gravity
p.runtime.h0 = 10000	# avg height
p.runtime.sphere_rotating_coriolis_omega = 7.292e-5	# coriolis effect
p.runtime.sphere_radius = 6371220	# radius

p.runtime.bench_id = -1
p.runtime.benchmark_name = "gaussian_bumps2"



####################################
# REXI
####################################

p.runtime.rexi_method = 'terry'

reference_timestep_size = 50
timestep_sizes = [50, 100, 200, 400, 800]



if __name__ == "__main__":

	####################################
	# RK 4 reference
	####################################

	if True:
		p.runtime.timestepping_method = 'l_erk'
		p.cluster.pm_time_cores_per_mpi_rank = 1

		for p.runtime.timestepping_order in [4]:
			for p.runtime.timestep_size in [50]:
			#for p.runtime.timestep_size in timestep_sizes:
				p.gen_script('script_ref'+p.runtime.getUniqueID(p.compile), 'run.sh')


	####################################
	# RK
	####################################

	if True:
		p.runtime.timestepping_method = 'l_erk'
		p.cluster.pm_time_cores_per_mpi_rank = 1

		for p.runtime.timestepping_order in [2]:
			#for p.runtime.timestep_size in [50]:
			for p.runtime.timestep_size in [50, 100, 200, 400, 800]:
				p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')




	####################################
	# Crank Nicolson
	####################################

	if True:
		p.runtime.timestepping_method = 'l_cn'
		p.cluster.pm_time_cores_per_mpi_rank = 1

		for p.runtime.timestepping_order in [2]:
			for p.runtime.timestep_size in [50, 100, 200, 400, 800]:
			#for p.runtime.timestep_size in timestep_sizes:
				p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')



	####################################
	# REXI short time steps
	####################################

	if True:
		p.runtime.timestepping_method = 'l_rexi'

		#for p.runtime.rexi_m in [128, 256, 512, 1024, 2048]:
		for p.runtime.rexi_m in [512, 1024]:
			for p.runtime.timestepping_order in [0]:
				for p.runtime.timestep_size in [800, 1600, 2*1600]:
				#for p.runtime.timestep_size in timestep_sizes:
					for p.cluster.par_time_cores in par_time_cores_list:
						p.gen_script('script'+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')



	####################################
	# REXI long time step
	####################################

	if True:
		p.runtime.timestepping_method = 'l_rexi'

		#for p.runtime.rexi_m in [1024, 4192]:
		for p.runtime.rexi_m in [4192]:
			for p.runtime.timestepping_order in [0]:
				for p.runtime.timestep_size in [129600]:
				#for p.runtime.timestep_size in timestep_sizes:
					for p.cluster.par_time_cores in par_time_cores_list:
						p.gen_script('script'+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')

