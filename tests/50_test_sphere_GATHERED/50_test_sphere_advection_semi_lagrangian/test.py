#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.unit_test="test_sphere_advection_semi_lagrangian"

jg.compile.plane_spectral_space="disable"
jg.compile.sphere_spectral_space="enable"
jg.runtime.benchmark_name = "adv_gauss_bump"

jg.runtime.semi_lagrangian_approximate_sphere_geometry = 1

jg.runtime.max_simulation_time = 60*60*24*12

jg.runtime.verbosity = 5

#params_compile_mode = ['release', 'debug']
params_compile_mode = ['release']
#params_compile_mode = ['debug']

params_runtime_mode_res_x = [64]
params_runtime_mode_res_y = [64]

params_advection_rotation_angles = [1.5708, 0, -0.7]
params_rotation_velocity = [0, 60*60*24*12]

params_runtime_ts_methods = [
        # [ method_id, order ]
        #["na_sl", 1],
        ["na_sl", 2],
        ["na_erk", 2],
#        ["na_erk", 4],
    ]


for (res_x, res_y) in product(params_runtime_mode_res_x, params_runtime_mode_res_y):
    jg.runtime.space_res_spectral = (res_x, res_y)

    # Iterate over time stepping methods and the order
    for ts_method in params_runtime_ts_methods:
        jg.runtime.timestepping_method = ts_method[0]
        jg.runtime.timestepping_order = ts_method[1]

        for vel in product(
            params_advection_rotation_angles,
            params_rotation_velocity,
        ):
            jg.runtime.advection_rotation_angle = vel[0]
            jg.runtime.advection_velocity = ",".join(str(x) for x in vel)

            jg.runtime.semi_lagrangian_max_iterations = 10

            if jg.runtime.timestepping_method == "na_sl":
                if jg.runtime.timestepping_order == 1:
                    params_runtime_timestep_sizes = [60*60*12]
                elif jg.runtime.timestepping_order == 2:
                    params_runtime_timestep_sizes = [60*60*24]
                else:
                    raise Exception("Unsupported order "+str(jg.runtime.timestepping_order))

            elif jg.runtime.timestepping_method == "na_erk":
                params_runtime_timestep_sizes = [60*60]
            else:
                raise Exception("Unknown time stepping method")

            for (
                jg.compile.mode,
                jg.runtime.timestep_size,
            ) in product(
                params_compile_mode,
                params_runtime_timestep_sizes,
            ):
                jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
