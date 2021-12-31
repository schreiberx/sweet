#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule_local.JobMule import *
jg = JobGeneration()

############################################################
# Settings that will stay the same for reference & test jobs
############################################################

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

jg.runtime.output_file_mode = 'bin'

# Verbosity mode
jg.runtime.verbosity = 2

# Mode and Physical resolution
jg.runtime.space_res_spectral = 64
jg.runtime.space_res_physical = None

# Benchmark
jg.runtime.benchmark_name = "galewsky"

# Compute error
jg.runtime.compute_error = 0

jg.runtime.f_sphere = 0

jg.runtime.viscosity = 0.0

jg.unique_id_filter = ['compile', 'parallelization']

timestep_size_min = 16
jg.runtime.max_simulation_time = timestep_size_min*1024
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

#####################################
# Reference Job: swe_sphere with ERK4
#####################################

jg.compile.program = 'swe_sphere'

ref_ts_size = 8
ref_ts_method = 'l_erk'
ref_ts_order = 4
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


#################################
# Test Jobs: libpfasst_swe_sphere
#################################

jg.compile.program = 'libpfasst_swe_sphere_mlsdc'
jg.compile.libpfasst = 'enable'

# LibPFASST runtime parameters
# set them all explicitly to make sure we know what's happening
jg.runtime.libpfasst_nlevels = 1
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

timestep_sizes = [timestep_size_min*(2.0**i) for i in range(0, 6)]

#
# Create job scripts
#

#for jg.runtime.libpfasst_nnodes in [3,5]:
for jg.runtime.libpfasst_niters in range(1,3):
    for jg.runtime.timestep_size in timestep_sizes:

        if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
            print("simtime: "+str(jg.runtime.max_simulation_time))
            print("timestep_size: "+str(jg.runtime.timestep_size))
            raise Exception("Invalid time step size (not remainder-less dividable)")
        
        jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 3)

        jg.gen_jobscript_directory()
