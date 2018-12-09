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
jg.compile.program = 'burgers'

jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'
jg.compile.sphere_spectral_space = 'disable'
jg.compile.sphere_spectral_dealiasing = 'disable'


# Enable quad math per default for CI REXI method
jg.compile.quadmath = 'enable'


# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = None
jg.runtime.space_res_physical = 128

#jg.runtime.benchmark_name = "gaussian_bumps_phi_vort_div"
jg.runtime.benchmark_name = "70"

# Viscosity
jg.runtime.viscosity = 0.01

# Simulation time
jg.runtime.max_simulation_time = 0.1

# Output data
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# Compute error
jg.runtime.compute_error = 0


# TODO: time step size
timestep_size_reference = 0.0001
timestep_sizes = [0.0001*(2.0**i) for i in range(0, 11)]

# Threading
jg.compile.threading = 'omp'

# Filter for unique job id
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

	ref_ts_method = 'l_direct'
	ref_ts_size = 0.0001
	ref_ts_order = 1

	ts_methods = [
			'l_irk',
			'l_erk',
		]


elif group == "l2":

	ts_order = 2

	ref_ts_method = 'l_direct'
	ref_ts_size = 0.0001
	ref_ts_order = 2

	ts_methods = [
			'l_irk',
			'l_erk',
		]


elif group == "ln1":

	ts_order = 1

	ref_ts_method = 'ln_cole_hopf'
	ref_ts_size = 0.0001
	ref_ts_order = 1

	ts_methods = [
			'ln_erk',
			'ln_imex',
		]


elif group == "ln2":

	ts_order = 2

	ref_ts_method = 'ln_cole_hopf'
	ref_ts_size = 0.0001
	ref_ts_order = 2

	ts_methods = [
			'ln_erk',
			'ln_imex',
		]

else:
	raise Exception("Group "+group+" not supported")


#elif group == "ln4":
#
#	ts_order = 4
#
#	ref_ts_method = 'ln_cole_hopf'
#	ref_ts_order = 0
#
#	ts_methods = [
#			'ln_erk',
#		]

jg.runtime.plane_domain_size = 1


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

		#if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
		#	print("simtime: "+str(jg.runtime.max_simulation_time))
		#	print("timestep_size: "+str(jg.runtime.timestep_size))
		#	raise Exception("Invalid time step size (not remainder-less dividable)")

		if 'rexi' in jg.runtime.timestepping_method:
			jg.runtime.rexi_method = 'ci'
		else:
			jg.runtime.rexi_method = None

		jg.gen_jobscript_directory()
