#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule.JobMule import *
jg = JobGeneration()

###################################################
# Compilation Settings for reference & tests jobs #
###################################################

jg.compile.program = 'libpfasst_swe_sphere_imex_sdc'

# enable libpfasst
jg.compile.libpfasst = 'enable'

# run simulation on sphere, not plane
jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'

# enable MPI
jg.compile.sweet_mpi = 'enable'

# other compilation settings
jg.compile.libsph = 'enable'
jg.compile.threading = 'off'
jg.compile.libfft = 'enable'
jg.compile.quadmath = 'enable'

###############################################
# Runtime Settings for reference & tests jobs #
###############################################

jg.runtime.output_file_mode = 'bin'

# Verbosity mode
jg.runtime.verbosity = 2

# Mode and Physical resolution
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_physical = None

# Benchmark
jg.runtime.benchmark_name = "galewsky"

# Compute error
jg.runtime.compute_errors = 0

jg.runtime.f_sphere = 0

jg.runtime.viscosity = 0.0

jg.unique_id_filter = ['compile', 'parallelization']

timestep_size_min = 16
jg.runtime.max_simulation_time = timestep_size_min*256
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# LibPFASST runtime parameters
# set them all explicitly to make sure we know what's happening
jg.runtime.libpfasst_nlevels = 1
jg.runtime.libpfasst_use_rk_stepper = 0

#################
# Reference Job #
#################

jg.runtime.libpfasst_nnodes = 5
jg.runtime.libpfasst_niters = 8
jg.runtime.libpfasst_nodes_type = 'SDC_GAUSS_LOBATTO'

ref_ts_size = 8

jg.runtime.timestep_size = ref_ts_size

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id


#############
# Test Jobs #
#############

timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

#
# Create job scripts
#

jg.runtime.libpfasst_nodes_type = 'SDC_GAUSS_LOBATTO'
for jg.runtime.libpfasst_nnodes in [3,5]:
    for jg.runtime.libpfasst_niters in range(1,5):
        for jg.runtime.timestep_size in timestep_sizes:

            if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
                print("simtime: "+str(jg.runtime.max_simulation_time))
                print("timestep_size: "+str(jg.runtime.timestep_size))
                raise Exception("Invalid time step size (not remainder-less dividable)")
            
            jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)

            jg.gen_jobscript_directory()

jg.runtime.libpfasst_nodes_type = 'SDC_GAUSS_LEGENDRE'
for jg.runtime.libpfasst_nnodes in [3,5]:
    for jg.runtime.libpfasst_niters in range(2,5):
        for jg.runtime.timestep_size in timestep_sizes:

            if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
                print("simtime: "+str(jg.runtime.max_simulation_time))
                print("timestep_size: "+str(jg.runtime.timestep_size))
                raise Exception("Invalid time step size (not remainder-less dividable)")
            
            jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * (jg.runtime.libpfasst_nnodes - 2))

            jg.gen_jobscript_directory()
