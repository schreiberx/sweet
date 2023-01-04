#! /usr/bin/env python3

import os
import sys
import math

from itertools import product


from mule_local.JobGeneration import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *
jg = JobGeneration()

verbose = False
#verbose = True

##################################################

jg.compile.compiler = 'intel'
jg.compile.sweet_mpi = 'enable'
jg.runtime.space_res_spectral = 256
jg.runtime.reuse_plans = 'load' # Load (if it exists)
jg.compile.threading = 'omp'

jg.parallelization.core_affinity = 'compact'
jg.parallelization.core_oversubscription = False

gen_reference_solution = False

#jg.runtime.max_simulation_time = 432000 #timestep_size_reference*(2**6)*10
jg.runtime.max_simulation_time = 720 #timestep_size_reference*(2**6)*10

params_timestep_sizes_explicit = [30]
#timestep_sizes_explicit = [10, 20, 30, 60, 120, 180]

params_timestep_sizes_implicit = [30]
#timestep_sizes_implicit = [60, 120, 180, 360, 480, 600, 720]

params_timestep_sizes_rexi = [30]
#timestep_sizes_rexi = [60, 120, 180, 240, 300, 360, 480, 600, 720]

# Parallelization
#params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_node]
#params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_node+1)]

params_ptime_num_cores_per_rank = [1]

unique_id_filter = []

unique_id_filter.append('runtime.disc_space')
unique_id_filter.append('runtime.reuse_plans')
unique_id_filter.append('runtime.simparams')
unique_id_filter.append('runtime.benchmark')
unique_id_filter.append('runtime.timestepping_order')
unique_id_filter.append('runtime.max_wallclock_time')

unique_id_filter.append('parallelization.mpi_ranks')
#unique_id_filter.append('parallelization.dims')
unique_id_filter.append('parallelization.cores_per_rank')
#unique_id_filter.append('parallelization.threads_per_rank')


jg.unique_id_filter = unique_id_filter


jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# No output
jg.runtime.output_filename = "-"

# REXI stuff
params_ci_N = [128]
#params_ci_N = [1]
params_ci_max_imag = [30.0]
params_ci_max_real = [10.0]

##################################################


#
# Force deactivating Turbo mode
#
jg.parallelization.force_turbo_off = True


def estimateWallclockTime(jg):
    # 30 minutes
    return 60*30

    """
    Return an estimated wallclock time
    """
    #
    # Reference wallclock time and corresponding time step size
    # e.g. for explicit RK2 integration scheme
    #
    # On Cheyenne with GNU compiler
    # OMP_NUM_THREADS=18
    # 247.378 seconds for ln_erk2 with dt=30 m=128 t=432000
    #
    ref_wallclock_seconds = 60*4
    ref_simtime = 432000
    ref_timestep_size = 60
    ref_mode_res = 128

    # Use this scaling for additional wallclock time
    safety_scaling = 10
    # 5 Min additionaly
    safety_add = 60*5

    wallclock_seconds = ref_wallclock_seconds

    # inv. linear with simulation time
    wallclock_seconds *= jg.runtime.max_simulation_time/ref_simtime

    # linear with time step size
    wallclock_seconds *= ref_timestep_size/jg.runtime.timestep_size

    # inverse quadratic with resolution
    wallclock_seconds *= pow(jg.runtime.space_res_spectral/ref_mode_res, 2.0)

    if jg.runtime.rexi_method != '':
        if jg.runtime.rexi_method != 'ci':
            raise Exception("TODO: Support other REXI methods")

        # Complex-valued
        wallclock_seconds *= 2.0

        # Number of REXI terms
        wallclock_seconds *= jg.runtime.rexi_ci_n

        # Parallelization in time
        wallclock_seconds /= jg.parallelization.pardims_dict['time'].num_ranks

    if wallclock_seconds <= 0:
        raise Exception("Estimated wallclock_seconds <= 0")

    wallclock_seconds *= safety_scaling
    wallclock_seconds += safety_add

    return wallclock_seconds

jg.compile.lapack = 'enable'
jg.compile.mkl = 'disable'

# Request dedicated compile script
jg.compilecommand_in_jobscript = False


#
# Run simulation on plane or sphere
#
jg.compile.program = 'swe_sphere'

jg.compile.plane_spectral_space = 'disable'
jg.compile.plane_spectral_dealiasing = 'disable'
jg.compile.sphere_spectral_space = 'enable'
jg.compile.sphere_spectral_dealiasing = 'enable'

jg.compile.benchmark_timings = 'enable'


#
# Activate Fortran source
#
jg.compile.fortran_source = 'enable'


# Verbosity mode
jg.runtime.verbosity = 0

#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = 128
jg.runtime.space_res_physical = -1

#
# Benchmark
#
jg.runtime.benchmark_name = "galewsky"

#
# Compute error
#
jg.runtime.compute_error = 0

#
# Preallocate the REXI matrices
#
jg.runtime.rexi_sphere_preallocation = 1

#
# Deactivate instability checks
#
jg.runtime.instability_checks = 0

jg.compile.rexi_thread_parallel_sum = 'disable'


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
#
jg.runtime.rexi_method = ''
jg.runtime.rexi_ci_n = 128
jg.runtime.rexi_ci_max_real = -999
jg.runtime.rexi_ci_max_imag = -999
jg.runtime.rexi_ci_sx = -1
jg.runtime.rexi_ci_sy = -1
jg.runtime.rexi_ci_mu = 0
jg.runtime.rexi_ci_primitive = 'circle'

#jg.runtime.rexi_beta_cutoff = 1e-16
#jg.runtime.rexi_beta_cutoff = 0

jg.runtime.viscosity = 0.0


timestep_size_reference = params_timestep_sizes_explicit[0]



#
# allow including this file
#
if __name__ == "__main__":

    # 2nd order nonlinear
    ts_methods = [
        ['ln_erk',        4,    4,    0],    # reference solution
        ['ln_erk',        2,    2,    0],

        #['lg_irk_lc_n_erk_ver0',    2,    2,    0],
        ['lg_irk_lc_n_erk_ver1',    2,    2,    0],
        #['l_irk_n_erk_ver0',    2,    2,    0],
        #['l_irk_n_erk_ver1',    2,    2,    0],

        #['lg_rexi_lc_n_erk_ver0',    2,    2,    0],
        #['lg_rexi_lc_n_erk_ver1',    2,    2,    0],
        #['l_rexi_n_erk_ver0',    2,    2,    0],
        #['l_rexi_n_erk_ver1',    2,    2,    0],

        #['lg_rexi_lc_n_etdrk',    2,    2,    0],
        #['l_rexi_n_etdrk',    2,    2,    0],
    ]

    #
    # Reference solution
    #
    if gen_reference_solution:
        tsm = ts_methods[0]
        jg.runtime.timestep_size  = timestep_size_reference

        jg.runtime.timestepping_method = tsm[0]
        jg.runtime.timestepping_order = tsm[1]
        jg.runtime.timestepping_order2 = tsm[2]

        # Update TIME parallelization
        ptime = JobParallelizationDimOptions('time')
        ptime.num_cores_per_rank = 1
        ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
        ptime.num_ranks = 1

        pspace = JobParallelizationDimOptions('space')
        pspace.num_cores_per_rank = 1
        pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]
        pspace.num_ranks = 1

        # Setup parallelization
        jg.setup_parallelization([pspace, ptime])

        if verbose:
            pspace.print()
            ptime.print()
            jg.parallelization.print()

        if len(tsm) > 4:
            s = tsm[4]
            jg.load_from_dict(tsm[4])

        jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

        jg.gen_jobscript_directory('job_benchref_'+jg.getUniqueID(unique_id_filter))



    #
    # Create job scripts
    #
    for tsm in ts_methods[1:]:

        jg.runtime.timestepping_method = tsm[0]
        jg.runtime.timestepping_order = tsm[1]
        jg.runtime.timestepping_order2 = tsm[2]

        if len(tsm) > 4:
            s = tsm[4]
            jg.runtime.load_from_dict(tsm[4])

        tsm_name = tsm[0]
        if 'ln_erk' in tsm_name:
            timestep_sizes = params_timestep_sizes_explicit
        elif 'l_erk' in tsm_name or 'lg_erk' in tsm_name:
            timestep_sizes = params_timestep_sizes_explicit
        elif 'l_irk' in tsm_name or 'lg_irk' in tsm_name:
            timestep_sizes = params_timestep_sizes_implicit
        elif 'l_rexi' in tsm_name or 'lg_rexi' in tsm_name:
            timestep_sizes = params_timestep_sizes_rexi
        else:
            print("Unable to identify time stepping method "+tsm_name)
            sys.exit(1)

        for pspace_num_cores_per_rank, pspace_num_threads_per_rank, jg.runtime.timestep_size in product(params_pspace_num_cores_per_rank, params_pspace_num_threads_per_rank, timestep_sizes):
            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = pspace_num_cores_per_rank
            pspace.num_threads_per_rank = pspace_num_threads_per_rank
            pspace.num_ranks = 1
            pspace.setup()


            if not '_rexi' in jg.runtime.timestepping_method:
                jg.runtime.rexi_method = ''

                # Update TIME parallelization
                ptime = JobParallelizationDimOptions('time')
                ptime.num_cores_per_rank = 1
                ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
                ptime.num_ranks = 1
                ptime.setup()

                jg.setup_parallelization([pspace, ptime])

                if verbose:
                    pspace.print()
                    ptime.print()
                    jg.parallelization.print()

                jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

                jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())

            else:

                for ci_N, ci_max_imag, ci_max_real in product(params_ci_N, params_ci_max_imag, params_ci_max_real):
                    jg.runtime.load_from_dict({
                        'rexi_method': 'ci',
                        'ci_n':ci_N,
                        'ci_max_real':ci_max_real,
                        'ci_max_imag':ci_max_imag,
                        'half_poles':0,
                        #'ci_gaussian_filter_scale':0.0,
                        #'ci_gaussian_filter_dt_norm':0.0,    # unit scaling for T128 resolution
                        #'ci_gaussian_filter_exp_N':0.0,
                    })

                    for time_ranks in params_ptime_num_cores_per_rank:

                            # Update TIME parallelization
                            ptime = JobParallelizationDimOptions('time')
                            ptime.num_cores_per_rank = 1
                            ptime.num_threads_per_rank = 1
                            ptime.num_ranks = time_ranks
                            ptime.setup()

                            jg.setup_parallelization([pspace, ptime])

                            if verbose:
                                pspace.print()
                                ptime.print()
                                jg.parallelization.print()

                            # Generate only scripts with max number of cores
                            jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

                            jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())


    if False:
        #
        # SHTNS plan generation scripts
        #
        jg.runtime.reuse_plans = 1    # search for awesome plans and store them

        #
        # Create dummy scripts to be used for SHTNS script generation
        #

        # No parallelization in time
        ptime = JobParallelizationDimOptions('time')
        ptime.num_cores_per_rank = 1
        ptime.num_threads_per_rank = 1
        ptime.num_ranks = 1
        ptime.setup()

        for tsm in ts_methods[1:2]:

            jg.runtime.timestepping_method = tsm[0]
            jg.runtime.timestepping_order = tsm[1]
            jg.runtime.timestepping_order2 = tsm[2]

            if not '_rexi' in jg.runtime.timestepping_method:
                jg.runtime.rexi_method = ''
            else:
                jg.runtime.rexi_method = 'ci'

            if len(tsm) > 4:
                s = tsm[4]
                jg.runtime.load_from_dict(tsm[4])


            for pspace_num_cores_per_rank, pspace_num_threads_per_rank, jg.runtime.timestep_size in product(params_pspace_num_cores_per_rank, params_pspace_num_threads_per_rank, timestep_sizes):
                pspace = JobParallelizationDimOptions('space')
                pspace.num_cores_per_rank = pspace_num_cores_per_rank
                pspace.num_threads_per_rank = pspace_num_threads_per_rank
                pspace.num_ranks = 1
                pspace.setup()

                jg.setup_parallelization([pspace, ptime])

                # Use 10 minutes per default to generate plans
                jg.parallelization.max_wallclock_seconds = 60*10

                # Set simtime to 0
                jg.runtime.max_simulation_time = 0

                # No output
                jg.runtime.output_timestep_size = -1
                jg.runtime.output_filename = "-"

                jg.gen_jobscript_directory('job_plan_'+jg.getUniqueID())



    # Write compile script
    jg.write_compilecommands("./compile_platform_"+jg.platforms.platform_id+".sh")


