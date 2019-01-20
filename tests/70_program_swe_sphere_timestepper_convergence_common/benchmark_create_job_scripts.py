#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule_local.JobMule import *
jg = JobGeneration()


#
# Run simulation on plane or sphere
#
jg.compile.program = 'swe_sphere'

jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'


# Enable quad math per default for CI REXI method
jg.compile.quadmath = 'enable'


# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_physical = None

jg.runtime.benchmark_name = "gaussian_bumps_phi_vort_div"

#
# Compute error
#
jg.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 1

#
# Threading accross all REXI terms
#
rexi_thread_par = True
if rexi_thread_par:
	# OMP parallel for over REXI terms
	jg.compile.threading = 'off'
	jg.compile.rexi_thread_parallel_sum = 'enable'
else:
	jg.compile.threading = 'omp'
	jg.compile.rexi_thread_parallel_sum = 'disable'



jg.runtime.f_sphere = 0

#jg.runtime.gravitation= 1
#jg.runtime.sphere_rotating_coriolis_omega = 1
#jg.runtime.h0 = 1
#jg.runtime.plane_domain_size = 1

jg.runtime.viscosity = 0.0


jg.unique_id_filter = ['compile', 'parallelization']

if len(sys.argv) <= 1:
	print("")
	print("Usage:")
	print("	"+sys.argv[0]+" [timestepping method]")
	print("")
	sys.exit(1)
	
group = sys.argv[1]





#
# allow including this file
#

if group == "l1":

	ts_order = 1

	ref_ts_method = 'l_erk'
	ref_ts_order = 4

	ts_methods = [
			#'l_erk_pvd',

			'l_erk',
			'lg_erk_lc_erk',

			'l_irk',
			'lg_irk_lc_erk',

			'l_rexi',
		]


elif group == "lg1":

	ts_order = 1

	ref_ts_method = 'lg_erk'
	ref_ts_order = 4

	ts_methods = [
			'lg_erk',
			'lg_irk',

			'lg_rexi',
		]

elif group == "l2":

	ts_order = 2

	ref_ts_method = 'l_erk'
	ref_ts_order = 4

	ts_methods = [

			'l_erk',
			'lg_erk_lc_erk',

			'l_cn',
			'lg_irk_lc_erk',

			'l_lf',
			'l_rexi',
		]

elif group == "lg2":

	ts_order = 2

	ref_ts_method = 'lg_erk'
	ref_ts_order = 4

	ts_methods = [
			'lg_erk',
			'lg_cn',

			'lg_rexi',
		]

elif group == "ln1":

	ts_order = 1

	ref_ts_method = 'ln_erk'
	ref_ts_order = 4

	ts_methods = [
			'ln_erk',

			'l_erk_n_erk',

			'lg_erk_lc_n_erk_ver0',
			'lg_erk_lc_n_erk_ver1',

			'l_irk_n_erk_ver0',
			'l_irk_n_erk_ver1',

			'lg_irk_lc_n_erk_ver0',
			'lg_irk_lc_n_erk_ver1',

			'l_rexi_n_erk_ver0',
			'l_rexi_n_erk_ver1',

			'lg_rexi_lc_n_erk_ver0',
			'lg_rexi_lc_n_erk_ver1',

			'l_rexi_n_etdrk',
		]


elif group == "ln2":

	ts_order = 2

	ref_ts_method = 'ln_erk'
	ref_ts_order = 4

	ts_methods = [
			'ln_erk',

			'l_erk_n_erk',

			'lg_erk_lc_n_erk_ver0',
			'lg_erk_lc_n_erk_ver1',

			'l_irk_n_erk_ver0',
			'l_irk_n_erk_ver1',

			'lg_irk_lc_n_erk_ver0',
			'lg_irk_lc_n_erk_ver1',

			'l_rexi_n_erk_ver0',
			'l_rexi_n_erk_ver1',

			'lg_rexi_lc_n_erk_ver0',
			'lg_rexi_lc_n_erk_ver1',

			'l_rexi_n_etdrk',
			'lg_rexi_lc_n_etdrk',

		]

else:
	raise Exception("Unknown group")


if ts_order == 1:
	ref_ts_size = 8
	timestep_size_min = 64
	timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

	jg.runtime.max_simulation_time = timestep_size_min*512
	jg.runtime.output_timestep_size = jg.runtime.max_simulation_time


elif ts_order == 2:
	#
	# A 2nd order accurate method already considerably reduces the errors
	# Therefore, we use larger time step sizes to increase the errors
	# to get errors larger than numerical precision
	#
	# We still want to have a very small time step size for the reference solution
	# This is in particular important for REXI comparisons with ln2-type tests
	#
	ref_ts_size = 8*2

	# Larger minimal time step size
	timestep_size_min = 64*4

	timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

	jg.runtime.max_simulation_time = timestep_size_min*512
	jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
	
else:
	raise Exception("Unsupported time integration order")

#
# CI REXI method
#
jg.runtime.rexi_method = 'ci'

# Use reduced number of REXI coefficients for convergence studies
if ts_order == 1:
	jg.runtime.rexi_ci_n = 16
	jg.runtime.rexi_ci_max_real = 1
	jg.runtime.rexi_ci_max_imag = 1

else:
	jg.runtime.rexi_ci_n = 32
	jg.runtime.rexi_ci_max_real = 2
	jg.runtime.rexi_ci_max_imag = 2


jg.runtime.rexi_ci_mu = 0
jg.runtime.rexi_ci_primitive = 'circle'



#
# Reference solution
#
jg.runtime.rexi_method = None
jg.runtime.timestepping_method = ref_ts_method
jg.runtime.timestepping_order = ref_ts_order
jg.runtime.timestepping_order2 = ref_ts_order
jg.runtime.timestep_size = ref_ts_size

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

		if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
			print("simtime: "+str(jg.runtime.max_simulation_time))
			print("timestep_size: "+str(jg.runtime.timestep_size))
			raise Exception("Invalid time step size (not remainder-less dividable)")

		if 'rexi' in jg.runtime.timestepping_method:
			jg.runtime.rexi_method = 'ci'
		else:
			jg.runtime.rexi_method = None

		jg.gen_jobscript_directory()
