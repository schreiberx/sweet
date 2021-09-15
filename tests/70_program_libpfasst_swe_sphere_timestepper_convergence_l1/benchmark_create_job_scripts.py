#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule_local.JobMule import *
jg = JobGeneration()

from mule_local.rexi.REXICoefficients import *
from mule_local.rexi.trexi.TREXI import *
from mule_local.rexi.cirexi.CIREXI import *
from mule_local.rexi.brexi.BREXI import *



#
# libpfasst_swe_sphere compilation options
#
jg.compile.program = 'libpfasst_swe_sphere'
jg.compile.libpfasst = 'enable'

# run simulation on sphere, not plane
jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'

# enable MPI
jg.compile.sweet_mpi = 'enable'

jg.compile.libsph = 'enable'
jg.compile.numa_block_allocator = 0
jg.compile.threading = 'off'

jg.compile.libfft = 'enable'

# Enable quad math per default for CI REXI method
jg.compile.quadmath = 'enable'


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
jg.runtime.compute_error = 0

jg.runtime.f_sphere = 0

#jg.runtime.gravitation= 1
#jg.runtime.sphere_rotating_coriolis_omega = 1
#jg.runtime.h0 = 1
#jg.runtime.plane_domain_size = 1

jg.runtime.viscosity = 0.0


jg.unique_id_filter = ['compile', 'parallelization']

# LibPFASST runtime parameters
# set them all explicitly to make sure we know what's happening
jg.runtime.libpfasst_nlevels = 1
#jg.runtime.libpfasst_niters = None
jg.runtime.libpfasst_nnodes = 5
jg.runtime.libpfasst_nsweeps_coarse = 1
jg.runtime.libpfasst_nodes_type = 'GAUSS_LOBATTO'
jg.runtime.libpfasst_coarsening_multiplier = 0.5
jg.runtime.libpfasst_use_rexi = 0
jg.runtime.libpfasst_implicit_coriolis_force = 0
jg.runtime.libpfasst_use_rk_stepper = 0


#####################################################
#####################################################
#####################################################

ref_nnodes = 5
ref_niters = 6
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
jg.runtime.timestep_size = ref_ts_size
jg.runtime.libpfasst_nnodes = ref_nnodes
jg.runtime.libpfasst_niters = ref_niters

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id



#
# Create job scripts
#

for jg.runtime.libpfasst_nnodes in [3,5]:
    for jg.runtime.libpfasst_niters in range(3,6):
        for jg.runtime.timestep_size in timestep_sizes:

            if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
                print("simtime: "+str(jg.runtime.max_simulation_time))
                print("timestep_size: "+str(jg.runtime.timestep_size))
                raise Exception("Invalid time step size (not remainder-less dividable)")

            jg.gen_jobscript_directory()
