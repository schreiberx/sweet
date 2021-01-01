#! /usr/bin/env python3

import os
import sys
import math

from itertools import product

# REXI
from mule_local.rexi.REXICoefficients import *
from mule_local.rexi.pcirexi.contour.BeanContour import BeanContour
from mule_local.rexi.pcirexi.contour.RectangleContour import RectangleContour
from mule_local.rexi.trexi.TREXI import *
from mule_local.rexi.cirexi.CIREXI import *
from mule_local.rexi.elrexi.ELREXI import *
from mule_local.rexi.brexi.BREXI import *
from mule_local.rexi.pcirexi.PCIREXI import *

# EFloat
efloat_mode = "float"


def get_rexi_benchmarks(jg):
    # Accumulator of all REXI methods
    # rexi_method['rexi_method'] = 'file'               # Choose REXI method which is typically 'file' for all file-based ones
    # rexi_method['rexi_files_coefficients'] = None     # List with approximations for different 'phi' functions
    rexi_methods = []

    #
    # CI REXI
    #
    if True:
        # REXI stuff
        def fun_params_ci_N(ci_max_real, ci_max_imag):
            if ci_max_imag >= 7:
                return 128
            else:
                return 32

        params_ci_max_imag = [30.0]
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

        rexi_method = {}

        # Choose REXI method which is typically 'file' for all file-based ones
        rexi_method['rexi_method'] = 'file'

        # List with approximations for different 'phi' functions
        rexi_method['rexi_files_coefficients'] = None

        for ci_max_imag, ci_max_real in product(params_ci_max_imag, params_ci_max_real):

            if params_ci_max_imag_scaling_relative_to_timestep_size != None:
                ci_max_imag *= (jg.runtime.timestep_size / params_ci_max_imag_scaling_relative_to_timestep_size)

            # "phi0"
            cirexi = CIREXI(efloat_mode=efloat_mode)
            coeffs_phi0 = cirexi.setup(function_name="phi0", N=fun_params_ci_N(ci_max_real, ci_max_imag),
                                       lambda_include_imag=ci_max_imag, lambda_max_real=ci_max_real).toFloat()

            # "phi1"
            cirexi = CIREXI(efloat_mode=efloat_mode)
            coeffs_phi1 = cirexi.setup(function_name="phi1", N=fun_params_ci_N(ci_max_real, ci_max_imag),
                                       lambda_include_imag=ci_max_imag, lambda_max_real=ci_max_real).toFloat()

            # "phi2"
            cirexi = CIREXI(efloat_mode=efloat_mode)
            coeffs_phi2 = cirexi.setup(function_name="phi2", N=fun_params_ci_N(ci_max_real, ci_max_imag),
                                       lambda_include_imag=ci_max_imag, lambda_max_real=ci_max_real).toFloat()

            rexi_method['rexi_files_coefficients'] = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

            # Add to list of REXI methods
            rexi_methods.append(rexi_method)

    #
    # EL-REXI
    #
    if True:
        max_imags = [30.0]
        rexi_method = {}

        # Choose REXI method which is typically 'file' for all file-based ones
        rexi_method['rexi_method'] = 'file'

        # List with approximations for different 'phi' functions
        rexi_method['rexi_files_coefficients'] = None

        for max_imag in max_imags:
            # "phi0"
            elrexi = ELREXI(efloat_mode=efloat_mode)
            coeffs_phi0 = elrexi.setup(function_name="phi0", N=max(64, int(75 * max_imag / 30)), lambda_max_real=10.5,
                                       lambda_max_imag=max_imag + 2.5).toFloat()

            # "phi1"
            elrexi = ELREXI(efloat_mode=efloat_mode)
            coeffs_phi1 = elrexi.setup(function_name="phi1", N=max(64, int(75 * max_imag / 30)), lambda_max_real=10.5,
                                       lambda_max_imag=max_imag + 2.5).toFloat()

            # "phi2"
            elrexi = ELREXI(efloat_mode=efloat_mode)
            coeffs_phi2 = elrexi.setup(function_name="phi2", N=max(64, int(75 * max_imag / 30)), lambda_max_real=10.5,
                                       lambda_max_imag=max_imag + 2.5).toFloat()

            rexi_method['rexi_files_coefficients'] = [coeffs_phi0, coeffs_phi1, coeffs_phi2]

            # Add to list of REXI methods
            rexi_methods.append(rexi_method)

    #
    # LR-REXI (Rectangle contour with Gauss-Legendre Quadrature)
    #
    if True:
        max_imags = [30.0]
        rexi_method = {}

        # Choose REXI method which is typically 'file' for all file-based ones
        rexi_method['rexi_method'] = 'file'

        # List with approximations for different 'phi' functions
        rexi_method['rexi_files_coefficients'] = None

        for max_imag in max_imags:
            # "phi0"
            pcirexi = PCIREXI()
            contour = RectangleContour(23, 2 * max_imag + 20, -1)
            i_s = InterpolationSettings(4, 2, 'equidistant')
            q_s = QuadratureSettings(overall_quadrature_points=128,
                                     quadrature_method="legendre",
                                     distribute_quad_points_based_on_arc_length=True)
            coeffs_phi0 = pcirexi.setup_pcirexi(contour=contour,
                                                interpolation_settings=i_s,
                                                quadrature_settings=q_s)
            coeffs_phi0.unique_id = "LR_REXI(width=23, height=" + str(2 * max_imag + 20) + ", center=-1, terms=128"
            rexi_method['rexi_files_coefficients'] = [coeffs_phi0]

            # Add to list of REXI methods
            rexi_methods.append(rexi_method)

    #
    # Bean-REXI
    #
    if True:
        max_imags = [30.0]
        rexi_method = {}

        # Choose REXI method which is typically 'file' for all file-based ones
        rexi_method['rexi_method'] = 'file'

        # List with approximations for different 'phi' functions
        rexi_method['rexi_files_coefficients'] = None

        for max_imag in max_imags:
            # "phi0"
            pcirexi = PCIREXI()
            contour = BeanContour(16, max_imag / 30 * 35, -2)
            i_s = InterpolationSettings(1, 0, 'direct')
            q_s = QuadratureSettings(overall_quadrature_points=max(64, int(75 * max_imag / 30)),
                                     quadrature_method="midpoint")
            coeffs_phi0 = pcirexi.setup_pcirexi(contour=contour,
                                                interpolation_settings=i_s,
                                                quadrature_settings=q_s)
            coeffs_phi0.unique_id = "Bean_REXI(width=16, height=" + str(
                max_imag / 30 * 35) + ", center=-2, terms=" + str(max(64, int(75 * max_imag / 30)))
            rexi_method['rexi_files_coefficients'] = [coeffs_phi0]

            # Add to list of REXI methods
            rexi_methods.append(rexi_method)

    return rexi_methods


if __name__ == "__main__":
    pass
