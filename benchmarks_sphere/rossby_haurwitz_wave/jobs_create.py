#! /usr/bin/env python3

import sys

from SWEET import *
p = SWEETJobGeneration()

p.cluster.setupTargetMachine('auto')


#
# Run simulation on plane or sphere
#
p.compile.program = 'swe_sphere'

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'disable'
p.compile.plane_spectral_dealiasing = 'disable'
p.compile.sphere_spectral_space = 'enable'
p.compile.sphere_spectral_dealiasing = 'enable'


#p.compile.compiler = 'intel'


#
# Use Intel MPI Compilers
#
#p.compile.compiler_c_exec = 'mpicc'
#p.compile.compiler_cpp_exec = 'mpicxx'
#p.compile.compiler_fortran_exec = 'mpif90'


#
# Activate Fortran source
#
p.compile.fortran_source = 'enable'


#
# MPI?
#
#p.compile.sweet_mpi = 'enable'


# Verbosity mode
p.runtime.verbosity = 2

#
# Mode and Physical resolution
#
p.runtime.space_res_spectral = 128
p.runtime.space_res_physical = -1

#
# Benchmark ID
# 4: Gaussian breaking dam
# 100: Galewski
#
p.runtime.benchmark_name = "rossby_haurwitz_wave"

#
# Compute error
#
p.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
#p.runtime.rexi_sphere_preallocation = 1

#
# Deactivate stability checks
#
p.stability_checks = 0

#
# Threading accross all REXI terms
#

p.compile.threading = 'off'


p.runtime.viscosity = 0.0


timestep_size_reference = 1

#timestep_sizes = [timestep_size_reference*(2**i) for i in range(2, 4)]
timestep_sizes = [1*(2**i) for i in range(0, 10)]

# 24 days
p.runtime.max_simulation_time = 24*24*60*60

# only one hour for convergence tests
p.runtime.max_simulation_time = 60*60

p.runtime.output_timestep_size = 60*60
#p.runtime.output_timestep_size = -1

p.runtime.rexi_extended_modes = 0

# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
groups = ['ln2']



p.cluster.par_space_cores = 16
p.cluster.pm_space_cores_per_mpi_rank = 16


#
# allow including this file
#
if __name__ == "__main__":

	groups = ['ln2']


	if len(sys.argv) > 1:
		groups = [sys.argv[1]]

	print("Groups: "+str(groups))

	for group in groups:

		# 2nd order nonlinear
		if group == 'ln2':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution

				###########
				#['l_irk_n_erk',		2,	2,	0],
				#['lg_irk_lc_n_erk_ver0',	2,	2,	0],
				['lg_irk_lc_n_erk_ver1',	2,	2,	0],
				['lg_erk_lc_n_erk',		2,	2,	0],
				#['l_cn_n_erk',			2,	2,	0],
				#['l_irk_n_erk',		2,	2,	0],
				#['l_rexi_n_erk',		2,	2,	0],
				#['l_rexi_n_etdrk',		2,	2,	0],


				#['ln_erk',		2,	2,	0],

				#['l_erk_n_erk',		2,	2,	0],

				#['lg_irk_lc_n_erk',	2,	2,	0],
				#['lg_rexi_lc_n_erk_ver0',	2,	2,	0],
				#['lg_rexi_lc_n_erk_ver1',	2,	2,	0],
				#['lg_rexi_lc_n_etdrk',	2,	2,	0],
			]

		# 4th order nonlinear
		if group == 'ln4':
			ts_methods = [
				['ln_erk',		4,	4,	0],	# reference solution
				['l_rexi_n_etdrk',	4,	4,	0],
				['ln_erk',		4,	4,	0],
			]




		#
		# add prefix string to group benchmarks
		#
		prefix_string_template = group



		#
		# Parallelization models
		#
		# Use 18 cores for each MPI task even if only 1 thread is used
		# This avoid any bandwidth-related issues
		#
		# Parallelization model (18 threads per rank)
#		p.cluster.pm_space_cores_per_mpi_rank = 18
#		p.cluster.pm_time_cores_per_mpi_rank = 1


		#
		# Reference solution
		#
		if True:
		#if False:
			tsm = ts_methods[0]
			p.runtime.timestep_size  = timestep_size_reference

			p.runtime.timestepping_method = tsm[0]
			p.runtime.timestepping_order = tsm[1]
			p.runtime.timestepping_order2 = tsm[2]
			p.runtime.rexi_use_direct_solution = tsm[3]

			p.cluster.par_time_cores = 1

			if len(tsm) > 4:
				s = tsm[4]
				p.load_from_dict(tsm[4])

			p.gen_script('script_'+prefix_string_template+'_ref'+p.runtime.getUniqueID(p.compile), 'run.sh')



		#
		# Create job scripts
		#
		for tsm in ts_methods[1:]:
			for p.runtime.timestep_size in timestep_sizes:
				p.runtime.timestepping_method = tsm[0]
				p.runtime.timestepping_order = tsm[1]
				p.runtime.timestepping_order2 = tsm[2]
				p.runtime.rexi_use_direct_solution = tsm[3]

				if len(tsm) > 4:
					s = tsm[4]
					p.runtime.load_from_dict(tsm[4])

				p.gen_script('script_'+prefix_string_template+p.runtime.getUniqueID(p.compile)+'_'+p.cluster.getUniqueID(), 'run.sh')



