#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule.JobMule import *
jg = JobGeneration()


#
# Run simulation on plane or sphere
#
jg.compile.program = 'swe_plane'

jg.compile.plane_spectral_space = 'enable'
jg.compile.plane_spectral_dealiasing = 'enable'
jg.compile.sphere_spectral_space = 'disable'
jg.compile.sphere_spectral_dealiasing = 'disable'

# Verbosity mode
jg.runtime.verbosity = 3

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = -1
jg.runtime.space_res_physical = 512

#
# Benchmark ID
# 1: Gaussian breaking dam
#
jg.runtime.benchmark_name = "rotated_steady_state"

#
# Compute error
#
# We have the real solution for these benchmarks
#
jg.runtime.compute_errors = 1

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 0

#
# Threading accross all REXI terms
#
rexi_thread_par = False
if rexi_thread_par:
    # OMP parallel for over REXI terms
    jg.compile.threading = 'off'
    jg.compile.rexi_thread_parallel_sum = 'enable'
else:
    jg.compile.threading = 'omp'
    jg.compile.rexi_thread_parallel_sum = 'disable'


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
# ==> use direct
jg.runtime.rexi_method = 'direct'
#jg.runtime.rexi_ci_n = 64
#jg.runtime.rexi_ci_sx = 50
#jg.runtime.rexi_ci_sy = 50
#jg.runtime.rexi_ci_mu = 0
#jg.runtime.rexi_ci_primitive = 'circle'
                
jg.runtime.gravitation= 1
jg.runtime.sphere_rotating_coriolis_omega = 1
jg.runtime.h0 = 100
jg.runtime.plane_domain_size = 1

jg.runtime.viscosity = 0.0


timestep_size_reference = 0.001
#timestep_sizes = [0.0001*(2.0**i) for i in range(0, 11)]

# Execute just a single time step
# Otherwise space interpolation errors are accumulating
jg.runtime.max_simulation_time = timestep_size_reference
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

phys_res_list = [16*(2**i) for i in range(0, 7)]

jg.unique_id_filter = ['compile', 'parallelization']



# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2']
groups = ['ln2space']


if len(sys.argv) > 1:
    groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
    # 1st order linear
    if group == 'l1':
        ref_ts_method = 'l_direct'
        ref_ts_order = 1
        ref_ts_size = jg.max_simulation_time

        ts_methods = [
            ['l_erk',    1,    0],
            ['l_irk',    1,    0],
            ['l_rexi',    1,    0],
        ]

    # 2nd order linear
    if group == 'l2':
        ref_ts_method = 'l_direct'
        ref_ts_order = 2
        ref_ts_size = jg.max_simulation_time

        ts_methods = [
            ['l_erk',    2,    0],
            ['l_cn',    2,    0],
            ['l_rexi',    2,    0],
        ]

    # 1st order nonlinear
    if group == 'ln1':
        ref_ts_method = 'ln_erk'
        ref_ts_order = 4
        ref_ts_size = timestep_size_reference

        ts_methods = [
            ['ln_erk',        1,    1],
            ['l_erk_n_erk',        1,    1],
            ['l_irk_n_erk',        1,    1],
            ['l_rexi_n_erk',    1,    1],
        ]

    # 1st order nonlinear
    if group == 'ln1test':
        ref_ts_method = 'ln_erk'
        ref_ts_order = 4
        ref_ts_size = timestep_size_reference

        ts_methods = [
            ['l_erk_n_erk',        1,    1],
            ['l_irk_n_erk',        1,    1],
            ['ln_erk',        1,    1],
        ]

    # 2nd order nonlinear
    if group == 'ln2':
        ref_ts_method = 'ln_erk'
        ref_ts_order = 4
        ref_ts_size = timestep_size_reference

        ts_methods = [
            ['l_cn_n_erk',        2,    2],
            ['l_erk_n_erk',        2,    2],
            ['l_irk_n_erk',        2,    2],
            ['ln_erk',        2,    2],
            ['l_rexi_n_erk',    2,    2],
        ]

    # 2nd order nonlinear non-fully-spectral
    if group == 'ln2space':
        ref_ts_method = 'ln_erk'
        ref_ts_order = 4
        ref_ts_size = timestep_size_reference

        ts_methods = [
            ['ln_erk',        2,    2],    # FD- C-grid
            ['l_cn_na_sl_nd_settls', 2,    2],    # SI-SL-SP
            ['l_rexi_na_sl_nd_settls',    2,    2], #SL-EXP-SETTLS
            ['l_rexi_na_sl_nd_etdrk',    2,    2], #SL-EXP-ETDRK
#            ['l_rexi_n_erk',    2,    2],
        ]

    #
    # OVERRIDE TS methods
    #
    if len(sys.argv) > 4:
        ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4])]]


    #
    # add prefix string to group benchmarks
    #
    prefix_string_template = group


    #
    # Reference solution
    # (Not required, since we compare it with the real errors with --compute-errors=1)
    #if True:
    if False:
        print("Reference")
        tsm = ts_methods[0]
    
        jg.timestep_size = timestep_size_reference
        jg.runtime.timestepping_method = tsm[0]
        jg.runtime.timestepping_order = tsm[1]
        jg.runtime.timestepping_order2 = tsm[2]
        jg.space_res_physical = 512

        if len(tsm) > 4:
            s = tsm[4]
            jg.runtime.load_from_dict(tsm[4])

        jg.gen_jobscript_directory()


    for tsm in ts_methods:

        if group == 'ln2space' and 'ln_erk' in tsm[0]:
            jg.runtime.space_grid_use_c_staggering = 1
            jg.runtime.space_use_spectral_basis_diffs = 0

            #jg.compile.plane_spectral_space = 'disable'
            jg.compile.plane_spectral_dealiasing = 'disable'
            jg.compile.libfft = 'enable'
        else:
            jg.runtime.space_grid_use_c_staggering = 0
            jg.runtime.space_use_spectral_basis_diffs = 1

            jg.compile.plane_spectral_space = 'enable'
            jg.compile.plane_spectral_dealiasing = 'enable'

        for phys_res in phys_res_list:

            jg.runtime.timestep_size = timestep_size_reference
            jg.runtime.timestepping_method = tsm[0]
            jg.runtime.timestepping_order = tsm[1]
            jg.runtime.timestepping_order2 = tsm[2]
            jg.runtime.space_res_physical = phys_res


            if len(tsm) > 4:
                s = tsm[4]
                jg.runtime.load_from_dict(tsm[4])

            jg.gen_jobscript_directory()

