#! /usr/bin/env python3

import os
import sys
import math

from itertools import product

# REXI
from mule.rexi.REXICoefficients import *
from mule.rexi.pcirexi.BeanREXI import BeanREXI
from mule.rexi.pcirexi.LRREXI import LRREXI
from mule.rexi.trexi.TREXI import *
from mule.rexi.cirexi.CIREXI import *
from mule.rexi.elrexi.ELREXI import *
from mule.rexi.brexi.BREXI import *

# EFloat
efloat_mode = "float"




def new_rexi_method():
    rexi_method = {}

    # Choose REXI method which is typically 'file' for all file-based ones
    rexi_method['rexi_method'] = 'file'

    # List with approximations for different 'phi' functions
    rexi_method['rexi_files_coefficients'] = None

    return rexi_method


def get_rexi_benchmarks(jg):
    # Accumulator of all REXI methods
    # rexi_method['rexi_method'] = 'file'               # Choose REXI method which is typically 'file' for all file-based ones
    # rexi_method['rexi_files_coefficients'] = None     # List with approximations for different 'phi' functions
    rexi_methods = []

    #for N in [64, 96, 128]:
    for N in [64]:
        #
        # CI REXI
        #
        if False:


            # REXI stuff
            def fun_params_ci_N(ci_max_real, ci_max_imag):
                raise Exception("XXX")
                return 64
                N = max(64, int(75 * ci_max_imag / 30))
                return N

            params_ci_max_imag = [30.0, 15.0]
            params_ci_max_real = [10.0]

            #
            # Scale the CI circle radius relative to this time step size
            # We do this simply to get a consistent time stepping method
            # Otherwise, CI would not behave consistently
            # Yes, that's ugly, but simply how it goes :-)
            #
            params_ci_max_imag_scaling_relative_to_timestep_size = 480
            # params_ci_max_imag_scaling_relative_to_timestep_size = None

            params_ci_min_imag = 5.0

            for ci_max_imag, ci_max_real in product(params_ci_max_imag, params_ci_max_real):
                rexi_method = new_rexi_method()

                #if params_ci_max_imag_scaling_relative_to_timestep_size != None:
                #    ci_max_imag *= (jg.runtime.timestep_size / params_ci_max_imag_scaling_relative_to_timestep_size)

                # "phi0"
                cirexi = CIREXI(efloat_mode=efloat_mode)
                coeffs_phi0 = cirexi.setup(
                        function_name="phi0",
                        N=N,
                        #N=fun_params_ci_N(ci_max_real, ci_max_imag),
                        #N=max(64, int(75 * ci_max_imag / 30)),
                        lambda_include_imag=ci_max_imag,
                        lambda_max_real=ci_max_real,
                        half_shifted = True
                    ).toFloat()

                # "phi1"
                cirexi = CIREXI(efloat_mode=efloat_mode)
                coeffs_phi1 = cirexi.setup(
                        function_name="phi1",
                        N=N,
                        #N=max(64, int(75 * ci_max_imag / 30)),
                        #N=fun_params_ci_N(ci_max_real, ci_max_imag),
                        lambda_include_imag=ci_max_imag,
                        lambda_max_real=ci_max_real,
                        half_shifted = True
                    ).toFloat()

                # "phi2"
                cirexi = CIREXI(efloat_mode=efloat_mode)
                coeffs_phi2 = cirexi.setup(
                        function_name="phi2",
                        N=N,
                        #N=max(64, int(75 * ci_max_imag / 30)),
                        #N=fun_params_ci_N(ci_max_real, ci_max_imag),
                        lambda_include_imag=ci_max_imag,
                        lambda_max_real=ci_max_real,
                        half_shifted = True
                    ).toFloat()

                rexi_method['rexi_files_coefficients'] = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                # Add to list of REXI methods
                rexi_methods.append(rexi_method)

        #
        # EL-REXI
        #
        if True:
            max_imags = [30.0]


            for normalize_rexi in [True, False]:
                for max_imag in max_imags:
                    rexi_method = new_rexi_method()

                    # "phi0"
                    elrexi = ELREXI(efloat_mode=efloat_mode)
                    coeffs_phi0 = elrexi.setup(
                            function_name="phi0",
                            #N=fun_params_ci_N(1, max_imag),
                            N=N,
                            #N=max(64, int(75 * max_imag / 30)),
                            lambda_max_real=10.5,
                            lambda_max_imag=max_imag + 2.5
                        ).toFloat()

                    # "phi1"
                    elrexi = ELREXI(efloat_mode=efloat_mode)
                    coeffs_phi1 = elrexi.setup(
                            function_name="phi1",
                            #N=fun_params_ci_N(1, max_imag),
                            N=N,
                            #N=max(64, int(75 * max_imag / 30)),
                            lambda_max_real=10.5,
                            lambda_max_imag=max_imag + 2.5
                        ).toFloat()

                    # "phi2"
                    elrexi = ELREXI(efloat_mode=efloat_mode)
                    coeffs_phi2 = elrexi.setup(
                            function_name="phi2",
                            #N=fun_params_ci_N(1, max_imag),
                            N=N,
                            #N=max(64, int(75 * max_imag / 30)),
                            lambda_max_real=10.5,
                            lambda_max_imag=max_imag + 2.5
                        ).toFloat()

                    if normalize_rexi:
                        coeffs_phi0.normalize_steady_state()
                        coeffs_phi1.normalize_steady_state()
                        coeffs_phi2.normalize_steady_state()

                    rexi_method['rexi_files_coefficients'] = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                    # Add to list of REXI methods
                    rexi_methods.append(rexi_method)

        #
        # LR-REXI (Rectangle contour with Gauss-Legendre Quadrature)
        #
        if False:
            max_imags = [30.0]

            for max_imag in max_imags:
                rexi_method = new_rexi_method()

                # "phi0"
                lrrexi = LRREXI(efloat_mode=efloat_mode)
                coeffs_phi0 = lrrexi.setup(
                        function_name="phi0",
                        width=23,
                        height=2 * max_imag + 20,
                        center=-1,
                        N=N,
                        #N=128
                    ).toFloat()

                # "phi1"
                lrrexi = LRREXI(efloat_mode=efloat_mode)
                coeffs_phi1 = lrrexi.setup(
                        function_name="phi1",
                        width=23,
                        height=2 * max_imag + 20,
                        center=-1,
                        N=N,
                        #N=128
                    ).toFloat()

                # "phi2"
                lrrexi = LRREXI(efloat_mode=efloat_mode)
                coeffs_phi2 = lrrexi.setup(
                        function_name="phi2",
                        width=23,
                        height=2 * max_imag + 20,
                        center=-1,
                        N=N,
                        #N=128
                    ).toFloat()

                rexi_method['rexi_files_coefficients'] = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                # Add to list of REXI methods
                rexi_methods.append(rexi_method)

        #
        # Bean-REXI
        #
        if False:
            max_imags = [30.0]

            for max_imag in max_imags:
                rexi_method = new_rexi_method()

                # "phi0"
                beanrexi = BeanREXI(efloat_mode=efloat_mode)
                coeffs_phi0 = beanrexi.setup(
                        function_name="phi0",
                        horizontal_radius=16,
                        vertical_radius=max_imag / 30 * 35,
                        center=-2,
                        #N=fun_params_ci_N(1, max_imag),
                        N=N,
                        #N=max(64, int(75 * max_imag / 30))
                    ).toFloat()

                # "phi1"
                beanrexi = BeanREXI(efloat_mode=efloat_mode)
                coeffs_phi1 = beanrexi.setup(
                        function_name="phi1",
                        horizontal_radius=16,
                        vertical_radius=max_imag / 30 * 35,
                        center=-2,
                        #N=fun_params_ci_N(1, max_imag),
                        N=N,
                        #N=max(64, int(75 * max_imag / 30))
                    ).toFloat()


                # "phi2"
                beanrexi = BeanREXI(efloat_mode=efloat_mode)
                coeffs_phi2 = beanrexi.setup(
                        function_name="phi2",
                        horizontal_radius=16,
                        vertical_radius=max_imag / 30 * 35,
                        center=-2,
                        #N=fun_params_ci_N(1, max_imag),
                        N=N,
                        #N=max(64, int(75 * max_imag / 30))
                    ).toFloat()

                rexi_method['rexi_files_coefficients'] = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

                # Add to list of REXI methods
                rexi_methods.append(rexi_method)

    return rexi_methods

