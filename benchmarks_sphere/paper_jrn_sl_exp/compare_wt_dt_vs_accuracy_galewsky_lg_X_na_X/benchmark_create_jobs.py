#! /usr/bin/env python3

import os
import sys
import math

from itertools import product

# REXI
from mule_local.rexi.REXICoefficients import *
from mule_local.rexi.trexi.TREXI import *
from mule_local.rexi.cirexi.CIREXI import *
from mule_local.rexi.brexi.BREXI import *

efloat_mode = "float"
#efloat_mode = "mpfloat"


from mule_local.JobGeneration import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *
jg = JobGeneration()

verbose = False
#verbose = True

##################################################

##################################################

jg.compile.mode = 'release'
if '_gnu' in os.getenv('MULE_PLATFORM_ID'):
    jg.compile.compiler = 'gnu'
else:
    jg.compile.compiler = 'intel'
jg.compile.sweet_mpi = 'enable'

jg.runtime.space_res_spectral = 128
#jg.runtime.reuse_plans = 2    # enforce using plans (todo, enforcing not yet implemented)!

jg.parallelization.core_oversubscription = False
jg.parallelization.core_affinity = 'compact'

jg.compile.threading = 'omp'
jg.compile.rexi_thread_parallel_sum = 'disable'

gen_reference_solution = True

jg.runtime.max_simulation_time = 60*60*24*5    # 5 days

#params_timestep_sizes_explicit = [30]
params_timestep_sizes_explicit = [15/8, 15/4, 15/2, 15, 30, 60, 120, 180, 360]

#params_timestep_sizes_implicit = [30]
params_timestep_sizes_implicit = [15/8, 15/4, 15/2, 15, 30, 60, 120, 180, 360, 480, 600, 720, 960]

#params_timestep_sizes_exp = [30]
params_timestep_sizes_exp = [15, 30, 60, 120, 180, 240, 300, 360, 480, 600, 720, 960]



# Parallelization
params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
#params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
params_ptime_num_cores_per_rank = [1]

unique_id_filter = []
#unique_id_filter.append('simparams')
unique_id_filter.append('compile')
unique_id_filter.append('disc_space')
unique_id_filter.append('timestep_order')
#unique_id_filter.append('timestep_size')
unique_id_filter.append('rexi_params')
unique_id_filter.append('benchmark')

jg.unique_id_filter = unique_id_filter


jg.runtime.output_timestep_size = jg.runtime.max_simulation_time

# No output
#jg.runtime.output_filename = "-"

# REXI stuff
def fun_params_ci_N(ci_max_real, ci_max_imag):
    if ci_max_imag >= 7:
        return 128
    else:
        return 32


#params_ci_N = [128]
#params_ci_N = [1]
params_ci_max_imag = [30.0]
params_ci_max_real = [10.0]

#
# Scale the CI circle radius relative to this time step size
# We do this simply to get a consistent time stepping method
# Otherwise, CI would not behave consistently
#
params_ci_max_imag_scaling_relative_to_timestep_size = 480
#params_ci_max_imag_scaling_relative_to_timestep_size = None

params_ci_min_imag = 5.0

##################################################


#
# Force deactivating Turbo mode
#
jg.parallelization.force_turbo_off = True


def estimateWallclockTime(jg):

    if jg.reference_job:
        return 2*24*60*60


    """
    Return an estimated wallclock time
    """


    if 'rexi' in jg.runtime.timestepping_method:
        return 12*60*60

    # This runtime is required for CN and explicit methods
    if jg.parallelization.num_ranks > 1:
        raise Exception("Shouldn't happen")

    # Give 1 hour for job at max
    return 1*60*60
    # Give 2h for non-REXI runs
    #return 2*60*60



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

    # quadratic with resolution
    wallclock_seconds *= pow(ref_mode_res/jg.runtime.space_res_spectral, 2.0)

    if jg.runtime.rexi_method != '':
        if jg.runtime.rexi_method != 'ci':
            raise Exception("TODO: Support other REXI methods")

        # Complex-valued
        wallclock_seconds *= 2.0

        # Number of REXI terms
        wallclock_seconds *= fun_params_ci_N(10, 10)

        # Parallelization in time
        wallclock_seconds /= jg.parallelization.pardims_dict['time'].num_ranks

    if wallclock_seconds <= 0:
        raise Exception("Estimated wallclock_seconds <= 0")

    wallclock_seconds *= safety_scaling
    wallclock_seconds += safety_add

    if wallclock_seconds > jg.platform_resources.max_wallclock_seconds:
        wallclock_seconds = jg.platform_resources.max_wallclock_seconds

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

jg.compile.quadmath = 'enable'


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

# Leave instability checks activated
# Don't activate them since they are pretty costly!!!
jg.runtime.instability_checks = 0


#
# REXI method
# N=64, SX,SY=50 and MU=0 with circle primitive provide good results
#
jg.runtime.rexi_method = ''

jg.runtime.viscosity = 0.0


timestep_size_reference = params_timestep_sizes_explicit[0]


jg.runtime.rexi_extended_modes = 0

# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
groups = ['ln2']





#
# allow including this file
#
if __name__ == "__main__":

    if len(sys.argv) > 1:
        groups = [sys.argv[1]]

    print("Groups: "+str(groups))

    for group in groups:
        # 1st order linear

        # 2nd order nonlinear
        if group == 'ln2':
            ts_methods = [
                ['l_na_erk',        4,    4,    0],    # reference solution

                ###########
                # Runge-Kutta
                ###########
                ['l_na_erk',        2,    2,    0],

                ###########
                # CN
                ###########
                #['lg_irk_lc_n_erk_ver0',    2,    2,    0],
                #['lg_irk_lc_n_erk_ver1',    2,    2,    0],

                #['l_irk_n_erk_ver0',    2,    2,    0],
                #['l_irk_n_erk_ver1',    2,    2,    0],

                ###########
                # REXI
                ###########
                #['lg_rexi_lc_n_erk_ver0',    2,    2,    0],
                #['lg_rexi_lc_n_erk_ver1',    2,    2,    0],

                #['l_rexi_n_erk_ver0',    2,    2,    0],
                #['l_rexi_n_erk_ver1',    2,    2,    0],

                ###########
                # ETDRK
                ###########
                #['lg_rexi_lc_n_etdrk',    2,    2,    0],
                #['l_rexi_n_etdrk',    2,    2,    0],

                ###########
                # SETTLS variants
                ###########
                #['l_rexi_n_erk_ver0',    2,    2,    0],
                ['l_irk_na_sl_settls',    2,    2,    0],
                ['lg_irk_na_sl_lc_settls',    2,    2,    0],
                ['lg_exp_na_sl_lc_settls',    2,    2,    0],

            ]


        # 4th order nonlinear
        if group == 'ln4':
            ts_methods = [
                ['ln_erk',        4,    4,    0],    # reference solution
                ['l_rexi_n_etdrk',    4,    4,    0],
                ['ln_erk',        4,    4,    0],
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

            jg.reference_job = True
            jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

            jg.gen_jobscript_directory('job_benchref_'+jg.getUniqueID())
            jg.reference_job = False

            jg.reference_job_unique_id = jg.job_unique_id



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

            exp_integrator = False

            tsm_name = tsm[0]
            if 'ln_erk' in tsm_name or 'l_na_erk' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_explicit
            elif 'l_erk' in tsm_name or 'lg_erk' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_explicit
            elif '_irk' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_implicit
            elif '_rexi' in tsm_name or '_exp' in tsm_name:
                params_timestep_sizes = params_timestep_sizes_exp
                exp_integrator = True
            else:
                raise Exception("Unable to identify time stepping method "+tsm_name)


            for pspace_num_cores_per_rank, pspace_num_threads_per_rank, jg.runtime.timestep_size in product(params_pspace_num_cores_per_rank, params_pspace_num_threads_per_rank, params_timestep_sizes):
                pspace = JobParallelizationDimOptions('space')
                pspace.num_cores_per_rank = pspace_num_cores_per_rank
                pspace.num_threads_per_rank = pspace_num_threads_per_rank
                pspace.num_ranks = 1
                pspace.setup()

                if not exp_integrator or 'lg_' in tsm_name:

                    # Always use direct REXI method if no parallel-in-time
                    jg.runtime.rexi_method = 'direct'

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

                    raise Exception("This branch shouldn't be taken, yet")

                    jg.runtime.rexi_method = 'file'
                    for ci_max_imag, ci_max_real in product(params_ci_max_imag, params_ci_max_real):

                        if params_ci_max_imag_scaling_relative_to_timestep_size != None:
                            ci_max_imag *= (jg.runtime.timestep_size/params_ci_max_imag_scaling_relative_to_timestep_size)

                        # "phi0"
                        cirexi = CIREXI(efloat_mode = efloat_mode)
                        coeffs_phi0 = cirexi.setup(function_name="phi0", N=fun_params_ci_N(ci_max_real, ci_max_imag), lambda_include_imag=ci_max_imag, lambda_max_real=ci_max_real).toFloat()

                        # "phi1"
                        cirexi = CIREXI(efloat_mode = efloat_mode)
                        coeffs_phi1 = cirexi.setup(function_name="phi1", N=fun_params_ci_N(ci_max_real, ci_max_imag), lambda_include_imag=ci_max_imag, lambda_max_real=ci_max_real).toFloat()

                        # "phi2"
                        cirexi = CIREXI(efloat_mode = efloat_mode)
                        coeffs_phi2 = cirexi.setup(function_name="phi2", N=fun_params_ci_N(ci_max_real, ci_max_imag), lambda_include_imag=ci_max_imag, lambda_max_real=ci_max_real).toFloat()

                        jg.runtime.rexi_files_coefficients = [coeffs_phi0, coeffs_phi1, coeffs_phi2]


                        # Update TIME parallelization
                        ptime = JobParallelizationDimOptions('time')
                        ptime.num_cores_per_rank = 1
                        ptime.num_threads_per_rank = 1
                        ptime.num_ranks = coeffs_phi0.len()
                        ptime.setup()

                        if jg.platform_resources.num_nodes == 1:
                            ptime.num_ranks = 1

                        jg.setup_parallelization([pspace, ptime])

                        if verbose:
                            pspace.print()
                            ptime.print()
                            jg.parallelization.print()

                        # Generate only scripts with max number of cores
                        jg.parallelization.max_wallclock_seconds = estimateWallclockTime(jg)

                        if int(jg.runtime.max_simulation_time / jg.runtime.timestep_size) * jg.runtime.timestep_size != jg.runtime.max_simulation_time:
                            raise Exception("Simtime "+str(jg.runtime.max_simulation_time)+" not dividable without remainder by "+str(jg.runtime.timestep_size))

                        jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())

    #
    # SHTNS plan generation scripts
    #
    #jg.runtime.reuse_plans = 1    # search for awesome plans and store them

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


        for pspace_num_cores_per_rank, pspace_num_threads_per_rank, jg.runtime.timestep_size in product(params_pspace_num_cores_per_rank, params_pspace_num_threads_per_rank, [params_timestep_sizes_explicit[0]]):
            pspace = JobParallelizationDimOptions('space')
            pspace.num_cores_per_rank = pspace_num_cores_per_rank
            pspace.num_threads_per_rank = pspace_num_threads_per_rank
            pspace.num_ranks = 1
            pspace.setup()

            jg.setup_parallelization([pspace, ptime])

            # Use 10 minutes per default to generate plans
            jg.parallelization.max_wallclock_seconds = 60*10

            # Set simtime to 0
            #jg.runtime.max_simulation_time = 0

            # No output
            jg.runtime.output_timestep_size = -1
            jg.runtime.output_filename = "-"

            jobdir = 'job_plan_'+jg.getUniqueID()
            jg.gen_jobscript_directory(jobdir)



    # Write compile script
    jg.write_compilecommands("./compile_platform_"+jg.platforms.platform_id+".sh")


