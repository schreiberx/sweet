#! /usr/bin/env python3

import os
import sys
import stat
import math

from mule.JobGeneration import *

#Create main compile/run options
jg = JobGeneration()

jg.compile.program = "swe_plane"

jg.runtime.verbosity = 2
jg.runtime.space_res_spectral = 128

jg.runtime.benchmark_name = 'benchmark_id_1'

jg.runtime.gravitation= 1
jg.runtime.sphere_rotating_coriolis_omega = 1
jg.runtime.h0 = 1
jg.runtime.plane_domain_size = 1

jg.runtime.rexi_method = 'direct'

jg.runtime.viscosity = 0.0

jg.runtime.max_simulation_time = 0.1
jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

jg.unique_id_filter = ['compile', 'runtime.benchmark', 'runtime.simparams']

timestep_size_reference = 0.0001
timestep_sizes = [0.0001*(2.0**i) for i in range(0, 11)]


# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
#groups = ['ln2test']

if len(sys.argv) > 1:
    groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
    #
    # Time integration method, order of 1st term, order of 2nd term
    #
    # 1st order linear
    if group == 'l1':
        ts_methods = [
            ['l_direct',    0,    0],    # reference solution
            ['l_erk',    1,    0],
            ['l_irk',    1,    0],
            ['l_rexi',    0,    0],
        ]

    # 2nd order linear
    if group == 'l2':
        ts_methods = [
            ['l_direct',    0,    0],    # reference solution
            ['l_erk',    2,    0],
            ['l_cn',    2,    0],
            ['l_rexi',    0,    0],
        ]

    #    ['lg_rexi_lc_erk_nt_sl_nd_erk',
    #    ['l_rexi_ns_sl_nd_erk',

    # 1st order nonlinear
    if group == 'ln1':
        ts_methods = [
            ['ln_erk',        4,    4],    # reference solution
            ['l_erk_n_erk',        1,    1],
            ['l_irk_n_erk',        1,    1],
            ['ln_erk',        1,    1],
            ['l_rexi_n_erk',    1,    1],
        ]

    # 2nd order nonlinear
    if group == 'ln2':
        ts_methods = [
            ['ln_erk',        4,    4],    # reference solution
            ['l_cn_n_erk',        2,    2],
            ['l_erk_n_erk',        2,    2],
            ['l_irk_n_erk',        2,    2],
            ['ln_erk',        2,    2],
            ['l_rexi_n_erk',    2,    2],
            ['l_rexi_ns_sl_nd_erk',    2,    2],
            ['lg_rexi_lc_erk_nt_sl_nd_erk',    2,    2],
        ]


    # 4th order nonlinear
    if group == 'ln4':
        ts_methods = [
            ['ln_erk',        4,    4],    # reference solution
            #['ln_etdrk',        4,    4],    # reference solution

            ['ln_etdrk',        4,    4],
            ['ln_erk',        4,    4],
        ]




    #
    # add prefix string to group benchmarks
    #
    prefix_string_template = group

    #
    # Reference solution
    #
    if True:
        print("")
        print("Reference solution")
        tsm = ts_methods[0]

        jg.runtime.timestepping_method = tsm[0]
        jg.runtime.timestepping_order = tsm[1]
        jg.runtime.timestepping_order2 = tsm[2]

        jg.runtime.timestep_size = timestep_size_reference

        # Tag this as a reference job
        jg.reference_job = True
        jg.gen_jobscript_directory()
        jg.reference_job = False

        # Use this one as the reference solution!
        jg.reference_job_unique_id = jg.job_unique_id


    for tsm in ts_methods[1:]:
        for jg.runtime.timestep_size in timestep_sizes:

            jg.runtime.timestepping_method = tsm[0]
            jg.runtime.timestepping_order = tsm[1]
            jg.runtime.timestepping_order2 = tsm[2]

            jg.gen_jobscript_directory()

