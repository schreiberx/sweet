#! /usr/bin/env python3

import os
import sys
import math

efloat_mode = "float"
#efloat_mode = "mpfloat"

from mule_local.JobGeneration import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *
jg = JobGeneration()

verbose = False
#verbose = True

##################################################
#   Software enviroment stuff
##################################################

jg.compile.mode = 'release'
if '_gnu' in os.getenv('MULE_PLATFORM_ID'):
    jg.compile.compiler = 'gnu'
else:
    jg.compile.compiler = 'intel'
jg.compile.sweet_mpi = 'disable'


jg.parallelization.core_oversubscription = False
jg.parallelization.core_affinity = 'compact'

jg.compile.threading = 'omp'
jg.compile.rexi_thread_parallel_sum = 'disable'

# Parallelization
params_pspace_num_cores_per_rank = [jg.platform_resources.num_cores_per_socket]
#params_pspace_num_threads_per_rank = [i for i in range(1, jg.platform_resources.num_cores_per_socket+1)]
params_pspace_num_threads_per_rank = [jg.platform_resources.num_cores_per_socket]
params_ptime_num_cores_per_rank = [1]

#
# Force deactivating Turbo mode
#
jg.parallelization.force_turbo_off = True

jg.compile.lapack = 'enable'
jg.compile.mkl = 'disable'

# Request dedicated compile script
jg.compilecommand_in_jobscript = False

##################################################
#   Basic compilation stuff
##################################################

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

jg.compile.fortran_source = 'enable'

##################################################
#   Basic simulation stuff
##################################################

jg.runtime.max_simulation_time = 60*60*24*365   # (in seconds) 1 year
#jg.runtime.max_simulation_time = 60*60*5    # 5 hours

jg.runtime.timestep_size = 120 # 1 minutes

jg.runtime.space_res_spectral = 128

jg.runtime.output_timestep_size = jg.runtime.max_simulation_time/730

# No output
#jg.runtime.output_filename = "output" #leave default
jg.runtime.output_file_mode = "csv_spec_evol"

# Verbosity mode
jg.runtime.verbosity = 1

#
# Benchmark
#
#jg.runtime.benchmark_name = "galewsky"
#jg.runtime.benchmark_name = "williamson2_linear"
jg.runtime.benchmark_name = "barotropic_vort_modes"

jg.runtime.timestepping_method = "lg_0_lc_n_erk_bv"
jg.runtime.timestepping_order = 3
jg.runtime.timestepping_order2 = 3

#
# Compute error
#
jg.runtime.compute_error = 0

# Leave instability checks activated
# Don't activate them since they are pretty costly!!!
jg.runtime.instability_checks = 0

jg.runtime.viscosity = 0.0

jg.runtime.rexi_method = ''

# Update TIME parallelization
ptime = JobParallelizationDimOptions('time')
ptime.num_cores_per_rank = 1
ptime.num_threads_per_rank = 1 #pspace.num_cores_per_rank
ptime.num_ranks = 1

pspace = JobParallelizationDimOptions('space')
pspace.num_cores_per_rank = 1
pspace.num_threads_per_rank = params_pspace_num_cores_per_rank[-1]/4
pspace.num_ranks = 1

jg.reference_job = False

# Setup parallelization
jg.setup_parallelization([pspace, ptime])

if verbose:
    pspace.print()
    ptime.print()
    jg.parallelization.print()

unique_id_filter = []
#unique_id_filter.append('runtime.disc_space')
unique_id_filter.append('runtime.simparams')
unique_id_filter.append('runtime.reuse_plans')
#unique_id_filter.append('runtime.benchmark') 
#unique_id_filter.append('runtime.timestepping') 
unique_id_filter.append('compile') 
unique_id_filter.append('parallelization') 
unique_id_filter.append('runtime.max_wallclock_time')
jg.unique_id_filter = unique_id_filter

#experiment paramters
sys.path.insert(1, '../')
import modes_experiment as mexp

#
# allow including this file
#
if __name__ == "__main__":


    basename = jg.runtime.benchmark_name

    more_samples = True
    if more_samples:
        alpha_min = 1.0
        alpha_max = 81.0
        alpha_samples = 80
    else:
        alpha_min = 1.0
        alpha_max = 20.0
        alpha_samples = 19

    n_list = [7, 3]
    m_list = [3, 1]
    n_list_back = [7, 3]
    m_list_back = [4, 0]
    experiment = mexp.modes_TC3(n_list, m_list, n_list_back, m_list_back, alpha_min, alpha_max, alpha_samples) 
    exp_filename = "mode_setup_n"+'-'.join(map(str, n_list))+"_m"+'-'.join(map(str, m_list))+".pckl"
    codes = experiment.codes
    
    experiment.save_file(exp_filename)
    
    
    #setup up mode initializations
    for mode_code in codes:
        jg.runtime.benchmark_name =  basename+"_"+mode_code
        jg.gen_jobscript_directory('job_bench_'+jg.getUniqueID())

    # Write compile script
    jg.write_compilecommands("./compile_platform_"+jg.platforms.platform_id+".sh")


