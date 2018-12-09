#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from SWEET import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = SWEETJobGeneration()
jg.compile.unit_test="test_plane_basic_advection"

jg.compile.plane_spectral_space="enable"
jg.runtime.benchmark_name = "gaussian_bump_advection"

jg.runtime.max_simulation_time = 10
jg.runtime.verbosity = 5

params_domain_size_scales = [1000000]

#params_compile_mode = ['release', 'debug']
params_compile_mode = ['release']

params_compile_plane_spectral_dealiasing = ['enable', 'disable']
#params_compile_plane_spectral_dealiasing = ['enable']

params_runtime_spectral_derivs = [0]

params_runtime_phys_res_x = [32]
params_runtime_phys_res_y = [32]

params_runtime_timestep_sizes = [0.1]

jg.runtime.h0 = 0

params_velocity_u = [20000]
params_velocity_v = [0]
params_advection_scheme = [2]
params_test_mode = [0, 1]
params_staggered_use_analytical_solution = [1]

params_runtime_ts_methods = [
		# [ method_id, order ]
		["nl_erk", 1],
		["nl_erk", 2],
		["nl_erk", 3],
		["nl_erk", 4],
	]


for (res_x, res_y) in product(params_runtime_phys_res_x, params_runtime_phys_res_y):
	jg.runtime.space_res_physical = (res_x, res_y)

	# Try out different variants of domain size
	for jg.runtime.plane_domain_size in product(params_domain_size_scales, params_domain_size_scales):

		# Iterate over time stepping methods and the order
		for ts_method in params_runtime_ts_methods:
			jg.runtime.timestepping_method = ts_method[0]
			jg.runtime.timestepping_order = ts_method[1]

			for (
				velocity_u,
				velocity_v,
				advection_scheme,
				test_mode,
				staggered_use_analytical_solution,
			) in product(
				params_velocity_u,
				params_velocity_v,
				params_advection_scheme,
				params_test_mode,
				params_staggered_use_analytical_solution,
			):
				jg.runtime.user_defined_parameters['vu'] = {'id': 'vu', 'value': velocity_u, 'option': '--velocity-u='}
				jg.runtime.user_defined_parameters['vv'] = {'id': 'vv', 'value': velocity_v, 'option': '--velocity-v='}
				jg.runtime.user_defined_parameters['as'] = {'id': 'as', 'value': advection_scheme, 'option': '--advection-scheme='}
				jg.runtime.user_defined_parameters['suas'] = {'id': 'suas', 'value': staggered_use_analytical_solution, 'option': '--staggered-use-analytical-solution='}
				jg.runtime.user_defined_parameters['tm'] = {'id': 'tm', 'value': test_mode, 'option': '--test-mode='}

				for (
					jg.compile.mode,
					jg.compile.plane_spectral_dealiasing,
					jg.runtime.space_use_spectral_basis_diffs,
					jg.runtime.timestep_size,
				) in product(
					params_compile_mode,
					params_compile_plane_spectral_dealiasing,
					params_runtime_spectral_derivs,
					params_runtime_timestep_sizes,
				):
					jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
	sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)
