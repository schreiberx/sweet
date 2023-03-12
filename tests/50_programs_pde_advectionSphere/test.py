#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="tests/pde_advectionSphere"

jg.runtime.benchmark_name = "williamson1b"

jg.runtime.max_simulation_time = 12*24*60*60
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

#jg.runtime.verbosity = 5

#jg.runtime.reuse_plans = "save"

jg.unique_id_filter = ['compile', 'parallelization', 'runtime.benchmark']


#params_compile_mode = ['debug']
params_compile_mode = ['release']

params_runtime_mode_res_x = [128]
params_runtime_mode_res_y = [128]

# rotation speeds
# 0: no rotation
# 20=simtime: one rotation
# 1: fast rotaitons
#params_rotation_velocity = [1, 20]
#params_advection_velocity_u = [0.1, -0.2]
#params_advection_velocity_v = [-0.1, 0.2]
#params_advection_velocity_u = [0.1, -0.2]
#params_advection_velocity_v = [0.2]

params_runtime_ts_methods = [
        # [ method_id, order, timestep size]
        ["na_sl", 1, 400],
        ["na_sl", 2, 400],
        ["na_erk", 2, 100],
        ["na_erk", 4, 100],
    ]


for (res_x, res_y) in product(params_runtime_mode_res_x, params_runtime_mode_res_y):
    jg.runtime.space_res_spectral = (res_x, res_y)

    # Iterate over time stepping methods and the order
    for ts_method in params_runtime_ts_methods:
        jg.runtime.timestepping_method = ts_method[0]
        jg.runtime.timestepping_order = ts_method[1]

        for jg.compile.mode in params_compile_mode:
            jg.runtime.timestep_size = ts_method[2]
            jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
