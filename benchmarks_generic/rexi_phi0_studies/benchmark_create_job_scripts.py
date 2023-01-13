#! /usr/bin/env python3

#from mule import *
from mule.JobMule import *
from mule.utils import exec_program

#
# REXI specific
#
from mule.rexi.REXICoefficients import *
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.brexi.BREXI import *

from itertools import product



jg = JobGeneration()

jg.compile.program = 'rexi_approximation'
jg.unique_id_filter = ['compile', 'parallelization']
jg.unique_id_filter +=  ['simparams', 'benchmark', 'timestep_order', 'timestep_size', 'disc_space', 'compile']

rexi_file_methods = ["trexi", "cirexi", "brexi"]

num_timesteps = 1
jg.runtime.timestep_size = 1.0

jg.runtime.user_defined_parameters['function_name'] = {
        'id': 'fun',
        'value' : 'phi0',
        'option' : '--function-name=',
    }

jg.runtime.user_defined_parameters['test_mode'] = {
        'id': 'tm',
        'value' : 0,
        'option' : '--test-mode=',
    }


lambda_params = numpy.array([numpy.power(2, float(n))*1.0j for n in range(-6, 8)])


#function_name_list = ["phi0", "phi1", "phi2", "phi3", "phi4", "ups1", "ups2", "ups3", "psi1", "psi2", "psi3"]
function_name_list = ["phi0"] #, "phi1", "phi2", "phi3", "phi4", "ups1", "ups2", "ups3"]


efloat_mode = "float"
#efloat_mode = "mpfloat"

for function_name in function_name_list:

    for lambda_val in lambda_params:

        jg.runtime.max_simulation_time = jg.runtime.timestep_size*num_timesteps

        jg.runtime.user_defined_parameters['lambda_real'] = {
                'id': 'lamre',
                'value' : lambda_val.real,
                'option' : '--lambda-real=',
            }

        jg.runtime.user_defined_parameters['lambda_imag'] = {
                'id': 'lamim',
                'value' : lambda_val.imag,
                'option' : '--lambda-imag=',
            }

        for rexi_file_method in rexi_file_methods:

            if rexi_file_method == "trexi":

                if function_name != "phi0":
                    continue

                trexi = TREXI(efloat_mode = efloat_mode)

                M_list = [64, 128, 256]
                h_list = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]

                for (M, h) in product(M_list, h_list):
                    coeffs = trexi.setup(M=M, h=h).toFloat()
                    jg.runtime.rexi_files_coefficients = [coeffs]
                    jg.gen_jobscript_directory()

                    if True:
                        # Validate with C-implementation of T-REXI method
                        jg.runtime.rexi_method = "terry"
                        jg.runtime.rexi_terry_m = M
                        jg.runtime.rexi_terry_h = h
                        jg.gen_jobscript_directory()

                        # Back to original version
                        jg.runtime.rexi_method = "file"


            elif rexi_file_method == "cirexi":

                cirexi = CIREXI(efloat_mode = efloat_mode)

                # CI-REXI: Number of quadrature poles
                N_list = [256, 512, 1024, 2048]

                # CI-REXI: Value on imaginary axis to be included
                lambda_include_imag_list = [15, 20, 40, 80]

                # CI-REXI: Maximum value of quadrature pole
                lambda_max_real_list = [1, 2, 4, 8, 16]

                for (N, lambda_include_imag, lambda_max_real) in product(N_list, lambda_include_imag_list, lambda_max_real_list):
                    coeffs = cirexi.setup(function_name=function_name, N=N, lambda_include_imag=lambda_include_imag, lambda_max_real=lambda_max_real).toFloat()
                    jg.runtime.rexi_files_coefficients = [coeffs]
                    jg.gen_jobscript_directory()

                    if True:
                        # Validate with C-implementation of CI-REXI method
                        jg.runtime.rexi_method = "ci"
                        jg.runtime.rexi_ci_n = N
                        jg.runtime.rexi_ci_max_real = lambda_max_real
                        jg.runtime.rexi_ci_max_imag = lambda_include_imag
                        jg.gen_jobscript_directory()

                        # Back to original version
                        jg.runtime.rexi_method = "file"


            elif rexi_file_method == "brexi":

                if function_name != "phi0":
                    continue

                brexi = BREXI(efloat_mode = efloat_mode)

                # CI-REXI: Number of quadrature poles
                N_list = [8, 10, 12, 14, 16]
                #N_list = [8, 10]

                quadrature_method_list = ["gauss", "radau", "chebyshev"]

                for (N, quadrature_method) in product(N_list, quadrature_method_list):
                    coeffs = brexi.setup(N=N, quadrature_method=quadrature_method).toFloat()
                    jg.runtime.rexi_files_coefficients = [coeffs]
                    jg.gen_jobscript_directory()

            else:
                raise Exception("Unknown method "+rexi_file_methods)


#exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
#if exitcode != 0:
#    sys.exit(exitcode)
#print("Benchmarks successfully finished")
#exec_program('mule.benchmark.cleanup_all', catch_output=False)

