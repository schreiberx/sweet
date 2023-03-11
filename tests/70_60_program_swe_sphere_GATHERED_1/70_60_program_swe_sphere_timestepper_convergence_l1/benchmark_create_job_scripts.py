#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule.JobMule import *
jg = JobGeneration()

from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *



#
# Run simulation on plane or sphere
#
jg.compile.program = 'programs/pde_sweSphere'

jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'


# Enable quad math per default for CI REXI method
#jg.compile.quadmath = 'enable'


# Verbosity mode
jg.runtime.verbosity = 2

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_physical = None

#jg.runtime.benchmark_name = "gaussian_bumps_test_cases"
jg.runtime.benchmark_name = "gaussian_bumps_pvd"

#
# Compute error
#
jg.runtime.compute_errors = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 1

#
# Threading accross all REXI terms
#
rexi_thread_par = True



jg.runtime.f_sphere = 0

#jg.runtime.gravitation= 1
#jg.runtime.sphere_rotating_coriolis_omega = 1
#jg.runtime.h0 = 1
#jg.runtime.plane_domain_size = 1

jg.runtime.viscosity = 0.0


jg.unique_id_filter = ['compile', 'parallelization']


#####################################################
#####################################################
#####################################################


ts_order = 1

ref_ts_method = 'l_erk'
ref_ts_order = 4

ts_methods = [
                'l_erk',
                'lg_erk_lc_erk',

                'l_irk',
                'lg_irk_lc_erk',

                'l_exp',
        ]

ref_ts_size = 8
timestep_size_min = 16
timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

jg.runtime.max_simulation_time = timestep_size_min*512

#####################################################
#####################################################
#####################################################

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time



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

        if '_exp' in jg.runtime.timestepping_method:

            if rexi_thread_par:
                # OMP parallel for over REXI terms
                jg.compile.threading = 'off'
                jg.compile.rexi_thread_parallel_sum = 'enable'
            else:
                jg.compile.threading = 'omp'
                jg.compile.rexi_thread_parallel_sum = 'disable'


            if 0:
                # CI REXI method in SWEET
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


            elif 0:
                # CI REXI via file input
                jg.runtime.rexi_method = 'file'

                cirexi = CIREXI()
                coeffs = cirexi.setup("phi0", N=32, R=2).toFloat()
                jg.runtime.rexi_files_coefficients = [coeffs]

            else:
                # B REXI via file
                jg.runtime.rexi_method = 'file'

                brexi = BREXI()
                coeffs = brexi.setup(N=8, quadrature_method='gauss_legendre').toFloat()
                jg.runtime.rexi_files_coefficients = [coeffs]

        else:
                jg.runtime.rexi_method = None

        jg.gen_jobscript_directory()
