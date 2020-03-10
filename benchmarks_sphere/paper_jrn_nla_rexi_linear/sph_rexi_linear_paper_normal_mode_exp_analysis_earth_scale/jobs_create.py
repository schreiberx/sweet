#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import sys

#from MULE import *
from mule_local.JobGeneration import JobGeneration
p = JobGeneration()

p.compile.compiler = 'intel'
p.compile.program = 'swe_sphere'
p.compile.fortran_source = 'enable'

p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'

p.compile.rexi_thread_parallel_sum = 'enable'
p.compile.threading = 'off'

p.runtime.space_res_spectral = 16
p.runtime.output_timestep_size = 0.01

p.runtime.gravitation= 1    # gravity
p.runtime.h0 = 100000    # avg height
p.runtime.sphere_rotating_coriolis_omega = 0.000145842    # coriolis effect
p.runtime.sphere_radius = 6371220    # radius

# 3: gaussian breaking dam
# 4: geostrophic balance test case
p.runtime.benchmark_name = "gaussian_bumps2"


p.runtime.max_simulation_time = 0.001 #math.inf

p.runtime.compute_error = 0


p.runtime.space_res_spectral = 16
p.runtime.normal_mode_analysis = 13

p.runtime.f_sphere = 1

p.runtime.rexi_sphere_preallocation = 0

default_timestep_size = 400
default_timesteps = 1 #default_simtime/default_timestep_size


for p.runtime.f_sphere in [1]:

    if p.runtime.f_sphere == -1:
        p.runtime.sphere_rotating_coriolis_omega = 0
        p.runtime.f_sphere = 0

    elif p.runtime.f_sphere == 0:
        # f-sphere
        p.runtime.sphere_rotating_coriolis_omega = 0.000072921    # \Omega coriolis effect

    else:
        #p.runtime.sphere_rotating_coriolis_omega = 0.000072921    # \Omega coriolis effect
        p.runtime.sphere_rotating_coriolis_omega = 2.0*p.runtime.f_sphere    # Constant f requires multiplication with 2.0


    ####################################
    # REXI dt=defaut_timestep_size
    ####################################
    if 1:

        p.runtime.rexi_method = "file"

        for p.runtime.rexi_extended_modes in [2]:
            p.runtime.timestepping_order = 0


            p.runtime.timestep_size = default_timestep_size
            p.runtime.max_simulation_time = default_timestep_size*default_timesteps
            p.runtime.max_timesteps_nr = default_timesteps

            for rexi_m in [2**i for i in range(4, 13)]:
                from mule_local.rexi.trexi.TREXI import *

                trexi = TREXI()
                coeffs = trexi.setup(M=rexi_m, h=0.2).toFloat()
                p.runtime.rexi_files_coefficients = [coeffs]

                p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')

    p.runtime.rexi_method = ''

    ####################################
    # RK1
    ####################################
    if 1:
        p.runtime.timestepping_method = 'l_erk'
        p.runtime.timestepping_order = 1

        p.runtime.timestep_size = default_timestep_size
        p.runtime.max_simulation_time = default_timestep_size*default_timesteps
        p.runtime.max_timesteps_nr = default_timesteps

        p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


    ####################################
    # RK2
    ####################################
    if 1:
        p.runtime.timestepping_method = 'l_erk'
        p.runtime.timestepping_order = 2

        p.runtime.timestep_size = default_timestep_size
        p.runtime.max_simulation_time = default_timestep_size*default_timesteps
        p.runtime.max_timesteps_nr = default_timesteps

        p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


    ####################################
    # RK4
    ####################################
    if 1:
        p.runtime.timestepping_method = 'l_erk'
        p.runtime.timestepping_order = 4

        p.runtime.timestep_size = default_timestep_size
        p.runtime.max_simulation_time = default_timestep_size*default_timesteps
        p.runtime.max_timesteps_nr = default_timesteps

        p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


    ####################################
    # IRK1
    ####################################
    if 0:
        p.runtime.timestepping_method = 'l_irk'
        p.runtime.timestepping_order = 1

        p.runtime.timestep_size = default_timestep_size
        p.runtime.max_simulation_time = default_timestep_size*default_timesteps
        p.runtime.max_timesteps_nr = default_timesteps

        p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


    ####################################
    # IRK2
    ####################################
    if 1:
        p.runtime.timestepping_method = 'l_cn'
        p.runtime.timestepping_order = 2

        p.runtime.timestep_size = default_timestep_size
        p.runtime.max_simulation_time = default_timestep_size*default_timesteps
        p.runtime.max_timesteps_nr = default_timesteps

        p.gen_script('script'+p.runtime.getUniqueID(p.compile), 'run.sh')


