#! /usr/bin/env python3

import sys

from SWEET import *
p = SWEETJobGeneration()

p.compile.compiler = 'intel'
p.compile.program = 'swe_sphere'
p.compile.fortran_source = 'enable'

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'

p.compile.rexi_thread_parallel_sum = 'enable'
p.compile.threading = 'off'




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
p.runtime.rexi_extended_modes = 0
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


####################################
# RK 4
####################################

if True:
	p.runtime.timestepping_method = 'l_erk'

	for p.runtime.timestepping_order in [4]:
		for p.runtime.timestep_size in [50]:
		#for p.runtime.timestep_size in timestep_sizes:
			p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


####################################
# RK 4
####################################

if True:
	p.runtime.timestepping_method = 'l_erk'

	for p.runtime.timestepping_order in [1, 2, 4]:
		#for p.runtime.timestep_size in [50]:
		for p.runtime.timestep_size in timestep_sizes:
			p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')



####################################
# Leapfrog
####################################

if True:
	p.runtime.timestepping_method = 'l_lf'

	for p.runtime.timestepping_order in [2]:
		#for p.runtime.timestep_size in [50]:
		for p.runtime.timestep_size in timestep_sizes:
			p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


####################################
# Crank Nicolson
####################################

if True:
	p.runtime.timestepping_method = 'l_cn'

	for p.runtime.timestepping_order in [2]:
		#for p.runtime.timestep_size in [50]:
		for p.runtime.timestep_size in timestep_sizes:
			p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')



####################################
# REXI
####################################

if True:
	p.runtime.timestepping_method = 'l_rexi'

	for p.runtime.rexi_m in [16, 32, 512]:
		for p.runtime.timestepping_order in [0]:
			for p.runtime.timestep_size in [800]:
			#for p.runtime.timestep_size in timestep_sizes:
				p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


####################################
# REXI short time steps
####################################

if True:
	p.runtime.timestepping_method = 'l_rexi'

	#for p.runtime.rexi_m in [16, 32, 64, 128, 256, 512]:
	for p.runtime.rexi_m in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]:
		for p.runtime.timestepping_order in [0]:
			for p.runtime.timestep_size in [800]:
			#for p.runtime.timestep_size in timestep_sizes:
				p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')



####################################
# REXI long time step
####################################

if True:
	p.runtime.timestepping_method = 'l_rexi'

	#for p.runtime.rexi_m in [1024, 4192]:
	for p.runtime.rexi_m in [128, 256, 512, 1024, 2048, 4192]:
		for p.runtime.timestepping_order in [0]:
			for p.runtime.timestep_size in [129600]:
			#for p.runtime.timestep_size in timestep_sizes:
				p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

