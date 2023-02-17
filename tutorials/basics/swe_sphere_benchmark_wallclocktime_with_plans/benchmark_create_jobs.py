#! /usr/bin/env python3

import os
import sys
import math

from itertools import product

from mule.JobGeneration import *

jg = JobGeneration()

#
# Compile options
#
jg.compile.program = 'swe_sphere'

#
# Runtime options
#
#jg.runtime.space_res_spectral = 256
jg.runtime.space_res_spectral = 64
#jg.runtime.max_simulation_time = 60*60*24*8    # 8 days
jg.runtime.max_simulation_time = 60*60*24*1    # 1 day
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time
jg.runtime.benchmark_name = "galewsky"

# Saves us some time in case of unstable simulations
jg.runtime.instability_checks = 1

# We require the plans to be recreated
# This allows better comparision of different runs since the plans could otherwise change between runs
jg.runtime.reuse_plans = "require_load"


#
# Which benchmarks to run
#

# Timestep sizes
#params_timestep_sizes_ = [15/8, 15/4, 15/2, 15, 30, 60, 120, 180, 360, 720]
params_timestep_sizes_ = [30, 60, 120, 180, 360, 720]

# Time stepping methods
# ["name", "order1", "order2"]
ts_methods = [
    ['ln_erk',        2,    2],
    ['lg_irk_lc_n_erk_ver1',    2,    2],
    ['l_irk_n_erk_ver1',    2,    2],
]



#
# Reference benchark
#

# Ts parameter
ref_ts_method = ['ln_erk',        4,    4]    # Used as reference solution

# Pick the smallest time step size for the reference time step size
timestep_size_reference = params_timestep_sizes_[0]



#
# Generate reference solution
#
jg.runtime.timestep_size  = timestep_size_reference

jg.runtime.timestepping_method = ref_ts_method[0]
jg.runtime.timestepping_order = ref_ts_method[1]
jg.runtime.timestepping_order2 = ref_ts_method[2]

# Set this to true to say that this is one of the reference jobs
jg.reference_job = True

jg.gen_jobscript_directory('job_benchref_'+jg.getUniqueID())
jg.reference_job = False

#
# Create job scripts
#

# We now reuse the unique job ID of the reference solution to tell the other jobs about their reference solution!
jg.reference_job_unique_id = jg.job_unique_id


#
# Create job scripts
#
for tsm in ts_methods:

    jg.runtime.timestepping_method = tsm[0]
    jg.runtime.timestepping_order = tsm[1]
    jg.runtime.timestepping_order2 = tsm[2]

    for jg.runtime.timestep_size in params_timestep_sizes_:
        jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())

#
# SHTNS plan generation scripts
#

#
# Search for plans and store them
# This is in particular important for running studies across several nodes
# since they rely on using the same transformation plans in order to have no
# load imbalances
#
jg.runtime.reuse_plans = "save"

#
# Create dummy scripts to be used for SHTNS script generation
#
for tsm in ts_methods:

    jg.runtime.timestepping_method = tsm[0]
    jg.runtime.timestepping_order = tsm[1]
    jg.runtime.timestepping_order2 = tsm[2]

    # Set simtime to 0
    jg.runtime.max_simulation_time = 0

    # No output
    jg.runtime.output_timestep_size = -1
    jg.runtime.output_filename = "-"

    jobdir = 'job_plan'+jg.getUniqueID()
    jg.gen_jobscript_directory(jobdir)

