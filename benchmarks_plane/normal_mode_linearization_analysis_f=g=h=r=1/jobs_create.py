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
#p.compile.plane_spectral_dealiasing = 'enable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'disable'
p.compile.sphere_spectral_dealiasing = 'disable'


p.compile.compiler = 'gnu'


# Verbosity mode
p.runtime.verbosity = 2

#
# Mode and Physical resolution
#
p.runtime.space_res_spectral = 32
p.runtime.space_res_physical = -1

#
# Benchmark ID
# 4: Gaussian breaking dam
#
#p.runtime.bench_id = 4

#
# Compute error
#
p.runtime.compute_error = 0


#
# Simulation parameters
#
p.runtime.sphere_rotating_coriolis_omega = 1
p.runtime.h0 = 1
p.runtime.gravitation= 1
p.runtime.a = 1


p.runtime.normal_mode_analysis = 1
p.runtime.timestep_size = 0.00001

default_timesteps = 1


####################################
# RK1
####################################
if 1:
	p.runtime.timestepping_method = 'l_erk'
	p.runtime.timestepping_order = 1

	p.runtime.simtime = p.runtime.timestep_size*default_timesteps
	p.runtime.max_timesteps = default_timesteps

	p.gen_script('script_'+p.runtime.getUniqueID(p.compile), 'run.sh')

