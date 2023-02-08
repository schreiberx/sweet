#! /usr/bin/env python3

from mule.JobGeneration import *
p = JobGeneration()

p.runtime.floating_point_output_digits = 16

p.compile.program = 'swe_sphere'
p.compile.mode = 'release'
#p.compile.mode = 'debug'

p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'

p.runtime.space_res_spectral = "66,66"
p.runtime.space_res_physical = "196,100"
p.runtime.instability_checks = 0
p.runtime.benchmark_name = 'galewsky'

for analytical in [0, 1]:

    #p.compile.threading = 'omp'
    p.compile.threading = 'off'

    p.unique_id_filter = ['runtime', 'compile', 'parallelization']

    #p.runtime.gravitation= 1
    #p.runtime.sphere_rotating_coriolis_omega = 1
    #p.runtime.h0 = 1
    #p.runtime.plane_domain_size = 1

    p.runtime.verbosity = 10

    if False:
            p.runtime.timestep_size = 60
            p.runtime.max_simulation_time = 60*10
            p.runtime.output_timestep_size = 60
    else:
            p.runtime.timestep_size = 10
            p.runtime.max_simulation_time = 60*60*24
            p.runtime.output_timestep_size = 60*60

    p.runtime.output_filename = "output_%s_t%020.8f.csv"

    if analytical:
        p.runtime.comma_separated_tags = "galewsky_analytical_geostrophic_setup"
        file_ext = "_analytical_1"
    else:
        p.runtime.comma_separated_tags = "galewsky_numerical_geostrophic_setup"
        file_ext = "_analytical_0"

    p.runtime.timestepping_method = 'ln_erk'
    p.runtime.timestepping_order = 1
    p.runtime.timestepping_order2 = 1

    p.gen_jobscript_directory('job_bench_sweet'+file_ext)

