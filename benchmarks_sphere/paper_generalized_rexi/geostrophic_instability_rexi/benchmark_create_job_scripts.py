#! /usr/bin/env python3

import sys

#from mule import *
from mule.JobMule import *
from mule.exec_program import *

#
# REXI specific
#
from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *

from itertools import product



jg = JobGeneration()

"""
Galwesky benchmark
"""
jg.compile.program = 'swe_sphere'

jg.compile.mode = 'release'
#p.compile.sweet_mpi = 'disable'


#
# Mode and Physical resolution
#
jg.runtime.space_res_spectral = 128
jg.runtime.space_res_physical = -1

jg.parallelization.core_oversubscription = False
jg.parallelization.core_affinity = 'compact'

jg.compile.threading = 'omp'
jg.compile.rexi_thread_parallel_sum = 'disable'

gen_reference_solution = False
jg.runtime.benchmark_name = "galewsky"

jg.runtime.max_simulation_time = 60*60*24*8    # 8 days

jg.runtime.output_timestep_size = 60*60  # Generate output every 1 hour
jg.runtime.output_file_mode = 'bin'

params_timestep_size_reference = 30.0

#params_timestep_sizes_explicit_ = [15*(2**i) for i in range(0, 4)]
#params_timestep_sizes_explicit_ = [60]

base_timestep_size = 128/jg.runtime.space_res_spectral*300.0
params_timestep_sizes_explicit_ = [base_timestep_size]

params_timestep_sizes_implicit_ = [15*(2**i) for i in range(2, 6)]
params_timestep_sizes_sl_ = [15*(2**i) for i in range(2, 6)]


"""
Setup filters to get nicer directory names
"""

jg.unique_id_filter = ['compile', 'parallelization']
jg.unique_id_filter +=  ['simparams', 'benchmark', 'timestep_order', 'timestep_size', 'disc_space', 'compile']

"""
REXI approximations
"""

# REXI methods
#rexi_file_methods = ["trexi", "cirexi", "brexi"]
rexi_file_methods = ["cirexi"]


# Number of time steps
num_timesteps = 10
jg.runtime.timestep_size = 1.0


#function_name_list = ["phi0", "phi1", "phi2", "phi3", "phi4", "ups1", "ups2", "ups3", "psi1", "psi2", "psi3"]
function_name_list = ["phi0"] #, "phi1", "phi2", "phi3", "phi4", "ups1", "ups2", "ups3"]


efloat_mode = "float"
#efloat_mode = "mpfloat"

#print(f"lambda_params: {lambda_params}")
#sys.exit(1)

jg.runtime.max_simulation_time = jg.runtime.timestep_size*num_timesteps

for rexi_file_method in rexi_file_methods:

    if rexi_file_method == "trexi":

        if function_name != "phi0":
            continue

        trexi = TREXI(efloat_mode = efloat_mode)

        M_list = [16, 32, 64, 128, 256, 512, 1024]
        h_list = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]

        for (M, h) in product(M_list, h_list):

            if len(function_name_list) != 1:
                raise Exception("Only phi0 allowed")

            if function_name_list[0] != "phi0":
                raise Exception("Only phi0 allowed")

            # Default poles
            coeffs = trexi.setup(M=M, h=h).toFloat()
            jg.runtime.rexi_files_coefficients = [coeffs]
            jg.gen_jobscript_directory()

            # Apply symmetric reduction
            coeffs.symmetric_reduction()
            jg.runtime.rexi_files_coefficients = [coeffs]
            jg.gen_jobscript_directory()


    elif rexi_file_method == "cirexi":

        cirexi = CIREXI(efloat_mode = efloat_mode)

        # CI-REXI: Number of quadrature poles
        #N_list = [32, 64, 128, 256, 512, 1024, 2048, 2*2048, 4*2048]
        N_list = [16]

        # CI-REXI: Value on imaginary axis to be included
        #lambda_include_imag_list = [5, 10, 15, 20, 40, 80, 120, 160, 200]
        lambda_include_imag_list = [10]

        # CI-REXI: Maximum value of quadrature pole along the real axis
        #lambda_max_real_list = [1, 2, 4, 8, 12, 16]
        lambda_max_real_list = [8]

        # beta filter
        beta_filter_list = [0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
        beta_filter_list = [0]

        for (N, lambda_include_imag, lambda_max_real, beta_filter) in product(N_list, lambda_include_imag_list, lambda_max_real_list, beta_filter_list):

            # Default poles
            coeffs_ = [cirexi.setup(function_name=function_name, N=N, lambda_include_imag=lambda_include_imag, lambda_max_real=lambda_max_real).toFloat() for function_name in function_name_list]
            jg.runtime.rexi_files_coefficients = coeffs_
            jg.gen_jobscript_directory()

            # Apply symmetric reduction
            for i in range(len(coeffs_)):
                coeffs_[i].symmetric_reduction()
            jg.runtime.rexi_files_coefficients = coeffs_
            jg.gen_jobscript_directory()

            if beta_filter != 0:
                # Apply beta filter after the symmetric reduction
                for i in range(len(coeffs_)):
                    coeffs_[i].beta_filter(beta_filter)
                jg.runtime.rexi_files_coefficients = coeffs_
                jg.gen_jobscript_directory()


    elif rexi_file_method == "brexi":

        if function_name != "phi0":
            continue

        brexi = BREXI(efloat_mode = efloat_mode)

        # CI-REXI: Number of quadrature poles
        N_list = [8, 10, 12, 14, 16, 32]
        #N_list = [8, 10]

        quadrature_method_list = ["gauss", "radau", "chebyshev"]

        for (N, quadrature_method) in product(N_list, quadrature_method_list):

            # Default poles
            coeffs = brexi.setup(N=N, quadrature_method=quadrature_method).toFloat()
            jg.runtime.rexi_files_coefficients = [coeffs]
            jg.gen_jobscript_directory()

            # Apply symmetric reduction
            coeffs.symmetric_reduction()
            jg.runtime.rexi_files_coefficients = [coeffs]
            jg.gen_jobscript_directory()


    else:
        raise Exception("Unknown method "+rexi_file_methods)


#exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
#if exitcode != 0:
#    sys.exit(exitcode)
#print("Benchmarks successfully finished")
#exec_program('mule.benchmark.cleanup_all', catch_output=False)

