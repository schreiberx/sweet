#! /usr/bin/env python3

import sys

from SWEET import *
p = JobGeneration()


#
# Cluster options
#
p.cluster.setupTargetMachine("cheyenne")
p.cluster.pm_space_cores_per_mpi_rank = 1
p.cluster.pm_time_cores_per_mpi_rank = 1
p.cluster.environment_vars = "export OMP_NUM_THREADS=1\n"
p.cluster.exec_prefix = " taskset -c 0 "



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
p.compile.sweet_mpi = 'disable'
#p.compile.sweet_mpi = 'enable'


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



default_ts_max_timesteps = 500

default_rexi_ci_n = 50
default_rexi_max_timesteps = default_ts_max_timesteps // default_rexi_ci_n

p.runtime.timestep_size = 1
p.runtime.max_simulation_time = 9999999999999999



for p.space_res_spectral in [512]:
	####################################
	# REXI
	####################################
	if True:
		p.runtime.timestepping_method = 'l_rexi'
		p.runtime.timestepping_order = 0
		p.runtime.rexi_method = 'ci'
		p.runtime.rexi_ci_n = default_rexi_ci_n
		p.runtime.rexi_ci_max_real = 10
		p.runtime.rexi_ci_max_imag = 10
		p.runtime.rexi_half_poles = 0
		p.runtime.sphere_extended_modes = 0
		p.runtime.max_timesteps_nr = default_rexi_max_timesteps

		p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

	p.runtime.rexi_method = ''
	p.runtime.max_timesteps_nr = default_ts_max_timesteps


	####################################
	# RKn
	####################################
	for i in [1, 2, 3, 4]:
		p.runtime.timestepping_method = 'l_erk'
		p.runtime.timestepping_order = i

		p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

	####################################
	# LF2
	####################################
	if True:
		p.runtime.timestepping_method = 'l_lf'
		p.runtime.timestepping_order = 2

		p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

	####################################
	# IRK1
	####################################
	for i in [1]:
		p.runtime.timestepping_method = 'l_irk'
		p.runtime.timestepping_order = i

		p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

	####################################
	# CN
	####################################
	for i in [2]:
		p.runtime.timestepping_method = 'l_cn'
		p.runtime.timestepping_order = i

		p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

