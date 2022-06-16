#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule_local.JobMule import *
jg = JobGeneration()

###################################################
# Compilation Settings for reference & tests jobs #
###################################################

jg.compile.program = 'libpfasst_swe_sphere_mlsdc'

# enable libpfasst
jg.compile.libpfasst = 'enable'

# run simulation on sphere, not plane
jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'

# enable MPI
jg.compile.sweet_mpi = 'enable'

jg.compile.libsph = 'enable'
jg.compile.threading = 'off'
jg.compile.libfft = 'enable'

# Enable quad math per default for CI REXI method
jg.compile.quadmath = 'enable'

jg.runtime.output_file_mode = 'bin'

# Verbosity mode
jg.runtime.verbosity = 2

# Mode and Physical resolution
jg.runtime.space_res_spectral = 256
jg.runtime.space_res_physical = None

# Benchmark
# jg.runtime.benchmark_name = "galewsky"
jg.runtime.benchmark_name = "gaussian_bump"

# Compute error
jg.runtime.compute_error = 0

jg.runtime.f_sphere = 0

jg.unique_id_filter = ['compile', 'parallelization']

jg.runtime.max_simulation_time = 86400
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# LibPFASST runtime parameters
# set them all explicitly to make sure we know what's happening
# jg.runtime.libpfasst_u2 = 1.0 * 1e4
jg.runtime.libpfasst_u2 = 1.0 * 1e5
jg.runtime.libpfasst_u4 = 0.0
jg.runtime.libpfasst_u6 = 0.0
jg.runtime.libpfasst_u8 = 0.0
jg.runtime.libpfasst_u_fields = "all"
jg.runtime.libpfasst_use_rk_stepper = 0

######################
# Reference Job: SDC #
######################

jg.runtime.libpfasst_nlevels = 1
jg.runtime.libpfasst_nnodes = 5
jg.runtime.libpfasst_niters = 8
jg.runtime.libpfasst_nodes_type = 'SDC_GAUSS_LOBATTO'

ref_ts_size = 90
jg.runtime.timestep_size = ref_ts_size

jg.reference_job = True
jg.gen_jobscript_directory()
jg.reference_job = False

# Use this one as the reference solution!
jg.reference_job_unique_id = jg.job_unique_id


#################################
# Test Jobs: libpfasst_swe_sphere
#################################

# SDC(2,2):

timestep_sizes = [15,30,60,120]

jg.runtime.libpfasst_nlevels = 1
jg.runtime.libpfasst_nnodes = 2
jg.runtime.libpfasst_niters = 2

for jg.runtime.timestep_size in timestep_sizes:
    if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
        print("simtime: "+str(jg.runtime.max_simulation_time))
        print("timestep_size: "+str(jg.runtime.timestep_size))
        raise Exception("Invalid time step size (not remainder-less dividable)")
    
    jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)
    jg.gen_jobscript_directory()

# SDC(3,4):

timestep_sizes = [60,120,240,400]

jg.runtime.libpfasst_nlevels = 1
jg.runtime.libpfasst_nnodes = 3
jg.runtime.libpfasst_niters = 4

for jg.runtime.timestep_size in timestep_sizes:
    if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
        print("simtime: "+str(jg.runtime.max_simulation_time))
        print("timestep_size: "+str(jg.runtime.timestep_size))
        raise Exception("Invalid time step size (not remainder-less dividable)")
    
    jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)
    jg.gen_jobscript_directory()

# SDC (5,8)

timestep_sizes = [120,240,400,800]

jg.runtime.libpfasst_nlevels = 1
jg.runtime.libpfasst_nnodes = 5
jg.runtime.libpfasst_niters = 8

for jg.runtime.timestep_size in timestep_sizes:
    if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
        print("simtime: "+str(jg.runtime.max_simulation_time))
        print("timestep_size: "+str(jg.runtime.timestep_size))
        raise Exception("Invalid time step size (not remainder-less dividable)")
    
    jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)
    jg.gen_jobscript_directory()

# MLSDC(3,2,2,1/8)

timestep_sizes = [15,30,60,120,240]

jg.runtime.libpfasst_nlevels = 2
jg.runtime.libpfasst_nnodes = 3
jg.runtime.libpfasst_niters = 2
jg.runtime.libpfasst_nsweeps_coarse = 2
jg.runtime.libpfasst_coarsening_multiplier = 0.125
jg.runtime.libpfasst_use_rk_stepper = 0

for jg.runtime.timestep_size in timestep_sizes:
    if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
        print("simtime: "+str(jg.runtime.max_simulation_time))
        print("timestep_size: "+str(jg.runtime.timestep_size))
        raise Exception("Invalid time step size (not remainder-less dividable)")
    
    jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)
    jg.gen_jobscript_directory()

# MLSDC(3,2,2,1/4)

timestep_sizes = [15,30,60,120,240]

jg.runtime.libpfasst_nlevels = 2
jg.runtime.libpfasst_nnodes = 3
jg.runtime.libpfasst_niters = 2
jg.runtime.libpfasst_nsweeps_coarse = 2
jg.runtime.libpfasst_coarsening_multiplier = 0.25
jg.runtime.libpfasst_use_rk_stepper = 0

for jg.runtime.timestep_size in timestep_sizes:
    if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
        print("simtime: "+str(jg.runtime.max_simulation_time))
        print("timestep_size: "+str(jg.runtime.timestep_size))
        raise Exception("Invalid time step size (not remainder-less dividable)")
    
    jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)
    jg.gen_jobscript_directory()

# MLSDC(3,2,2,1/2)

timestep_sizes = [30,60,120,240]

jg.runtime.libpfasst_nlevels = 2
jg.runtime.libpfasst_nnodes = 3
jg.runtime.libpfasst_niters = 2
jg.runtime.libpfasst_nsweeps_coarse = 2
jg.runtime.libpfasst_coarsening_multiplier = 0.5
jg.runtime.libpfasst_use_rk_stepper = 0

for jg.runtime.timestep_size in timestep_sizes:
    if jg.runtime.max_simulation_time % jg.runtime.timestep_size != 0:
        print("simtime: "+str(jg.runtime.max_simulation_time))
        print("timestep_size: "+str(jg.runtime.timestep_size))
        raise Exception("Invalid time step size (not remainder-less dividable)")
    
    jg.runtime.timestepping_order = min(jg.runtime.libpfasst_niters, 2 * jg.runtime.libpfasst_nnodes - 2)
    jg.gen_jobscript_directory()
