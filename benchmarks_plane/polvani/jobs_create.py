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

prefix=p.compile.program

p.compile.plane_or_sphere = 'sphere'
p.compile.plane_spectral_space = 'enable'
p.compile.plane_spectral_dealiasing = 'enable'
p.compile.sphere_spectral_space = 'disable'
p.compile.sphere_spectral_dealiasing = 'disable'


p.compile.compiler = 'gnu'


p.runtime.floating_point_output_digits = 12

#
# Use Intel MPI Compilers
#
#p.compile.compiler_c_exec = 'mpicc'
#p.compile.compiler_cpp_exec = 'mpicxx'
#p.compile.compiler_fortran_exec = 'mpif90'


#
# Activate Fortran source
#
#p.compile.fortran_source = 'enable'


#
# MPI?
#
#p.compile.sweet_mpi = 'enable'


# Verbosity mode
p.runtime.verbosity = 2

#
# Mode and Physical resolution
#
p.runtime.space_res_spectral = 200
p.runtime.space_res_physical = -1

#
# Benchmark ID
# 4: Gaussian breaking dam
# 100: Galewski
#
p.runtime.benchmark_name = 'polvani'

#
# Compute error
#
p.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
p.runtime.rexi_sphere_preallocation = 0

#
# Viscosity, hail to viscosity !!!
#
# Results in noisy output
#p.runtime.viscosity = 1e-10
#p.runtime.viscosity_order = 8

# Unstable!
#p.runtime.viscosity = 1e-8

# Stable
if p.runtime.space_res_spectral == 200:
	#p.runtime.viscosity = 1e-5
	p.runtime.viscosity = 2e-5
	p.runtime.viscosity_order = 8
else:
	p.runtime.viscosity = 1e-5*200.0/p.runtime.space_res_spectral
	p.runtime.viscosity_order = 8

#
# Deactivate stability checks
#
p.stability_checks = 1

#
# Threading accross all REXI terms

if True:
	p.compile.threading = 'omp'
	p.compile.rexi_thread_parallel_sum = 'disable'

else:
	#
	# WARNING: rexi_thread_par does not work yet!!!
	# MPI Ranks are clashing onthe same node with OpenMP Threads!
	#rexi_thread_par = True
	rexi_thread_par = False

	if rexi_thread_par:
		# OMP parallel for over REXI terms
		p.compile.threading = 'off'
		p.compile.rexi_thread_parallel_sum = 'enable'
	else:
		p.compile.threading = 'omp'
		p.compile.rexi_thread_parallel_sum = 'disable'


p.runtime.timestep_size = 0.005
p.runtime.timestep_size = 0.0025
p.runtime.max_simulation_time = 1000
p.runtime.output_timestep_size = 1000




#
# allow including this file
#
if __name__ == "__main__":
	p.runtime.timestepping_method = 'ln_erk'
	p.runtime.timestepping_order = 4

	if True:
	#for p.runtime.viscosity in [1.0*(1e-1**i) for i in range(0, 10)]:
		for [p.runtime.polvani_rossby, p.runtime.polvani_froude, M] in [
				[0.01, 0.04, 'A'],	# A
				[0.05, 0.05, 'B'],	# B
				[0.05, 0.075, 'C'],	# C
				[0.25, 0.05, 'D'],	# D
				[0.25, 0.20, 'E'],	# E
				[1.00, 0.05, 'F'],	# F
				[1.00, 0.30, 'G'],	# G
				[5.00, 0.05, 'H'],	# H
				[5.00, 0.30, 'I'],	# I
				[20.25, 0.30, 'J'],	# J
				[10.25, 0.10, 'K'],	# K
				[20.25, 0.05, 'L'],	# L
				[2.00, 0.10, 'M'],	# M
				[0.40, 0.10, 'N'],	# N
			]:
			p.gen_script('script_'+prefix+'_polvani_'+M+p.getUniqueID(), 'run.sh')

