#! /usr/bin/env python3

import os
import sys
import stat
import math

from SWEET import *
jg = SWEETJobGeneration()


#
# Run simulation on plane or sphere
#
jg.compile.program = 'swe_plane'

jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'
jg.compile.sphere_spectral_space = 'disable'
jg.compile.sphere_spectral_dealiasing = 'disable'

# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Physical resolution
#
jg.runtime.mode_res = -1
jg.runtime.phys_res = 128

#
# Benchmark ID
# 1: Gaussian breaking dam
#
jg.runtime.bench_id = 1

#
# Compute error
#
jg.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 0

#
# Threading accross all REXI terms
#
#rexi_thread_par = True
#if rexi_thread_par:
#	# OMP parallel for over REXI terms
#	jg.compile.threading = 'off'
#	jg.compile.rexi_thread_parallel_sum = 'enable'
#else:
#	jg.compile.threading = 'omp'
#	jg.compile.rexi_thread_parallel_sum = 'disable'

jg.compile.threading = 'omp'
jg.compile.rexi_thread_parallel_sum = 'disable'


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
#
#jg.runtime.rexi_method = 'ci'
jg.runtime.rexi_method = 'direct'
jg.runtime.rexi_ci_n = 64
jg.runtime.rexi_ci_max_real = 10
jg.runtime.rexi_ci_max_imag = 10
jg.runtime.rexi_ci_mu = 0
jg.runtime.rexi_ci_primitive = 'circle'

# Activate, since sometimes we need the numerical solution
jg.compile.quadmath = 'enable'


jg.runtime.g = 1
jg.runtime.f = 1
jg.runtime.h = 1
jg.runtime.domain_size = 1

jg.runtime.viscosity = 0.0


# This is the minimum for each ln2 time stepping methods to run below 10 minutes.
# This is important for Travis.
timestep_size_min = 0.0001

#
# This benchmark gets unstable for l_erk with dt=0.0016
#
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 5)]

# Don't use a smaller TS since convergence is not computable anymore

# Run for 500 time steps
jg.runtime.simtime = timestep_size_min*500


jg.runtime.output_timestep_size = jg.runtime.simtime



if len(sys.argv) <= 1:
	print("")
	print("Usage:")
	print("	"+sys.argv[0]+" [timestepping method]")
	print("")
	

group = sys.argv[1]


#
# allow including this file
#

if group == "l1":

	ts_order = 1
	jg.runtime.rexi_use_direct_solution = 1

	ref_ts_method = 'l_erk'
	ref_ts_order = 4
	ref_ts_size = timestep_size_min*0.5

	ts_methods = [
			'l_erk',
			'l_irk'
		]

elif group == "l2":

	ts_order = 2
	jg.runtime.rexi_use_direct_solution = 1

	ref_ts_method = 'l_erk'
	ref_ts_order = 4
	ref_ts_size = 0.0002

	ts_methods = [
			'l_erk',
			'l_cn',
		]


elif group == "ln1":

	ts_order = 1
	jg.runtime.rexi_use_direct_solution = 1

	ref_ts_method = 'ln_erk'
	ref_ts_order = 4
	ref_ts_size = timestep_size_min*0.5

	ts_methods = [
			'ln_erk',
			'l_erk_n_erk',

			'l_irk_n_erk',

			'l_rexi_n_erk',
			'l_rexi_n_etdrk',
		]


elif group == "ln2":

	ts_order = 2
	jg.runtime.rexi_use_direct_solution = 1

	ref_ts_method = 'ln_erk'
	ref_ts_order = 4
	ref_ts_size = 0.0002

	ts_methods = [
			'ln_erk',

			'l_cn_n_erk',
			'l_erk_n_erk',

			'l_rexi_n_erk',
			'l_rexi_n_etdrk',
		]

else:
	raise Exception("Unknown group")


#
# Reference solution
#
tsm = ts_methods[0]

jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = ref_ts_order
jg.runtime.timestepping_order2 = ref_ts_order
jg.runtime.timestep_size = ref_ts_size

if len(tsm) > 4:
	s = tsm[4]
	jg.runtime.load_from_dict(tsm[4])

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id



#
# Create job scripts
#

jg.runtime.timestepping_order = ts_order
jg.runtime.timestepping_order2 = ts_order

for tsm in ts_methods:
	for jg.runtime.timestep_size in timestep_sizes:
		jg.runtime.timestepping_method = tsm

		jg.gen_jobscript_directory()
