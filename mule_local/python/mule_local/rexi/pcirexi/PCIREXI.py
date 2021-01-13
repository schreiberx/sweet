#! /usr/bin/env python3
#
# Author: Raphael Schilling <raphael.schilling@tum.de>
#
import cmath

from mule_local.rexi.Functions import Functions
from mule_local.rexi.REXICoefficients import REXICoefficients
from mule_local.rexi.pcirexi.cauchy_integrator import CauchyIntegrator
from mule_local.rexi.pcirexi.contour.Contour import Contour
from mule_local.rexi.pcirexi.section.InterpolationSettings import InterpolationSettings
from mule_local.rexi.pcirexi.section.QuadratureSettings import QuadratureSettings
from mule_local.rexi.pcirexi.section.section_factory import SectionFactory
from mule_local.rexi.EFloat import *


class PCIREXI:
    contour: Contour = None
    interpolation_settings: InterpolationSettings = None
    quadrature_settings: QuadratureSettings = None

    #
    # Constructor
    # See setup(...) for documentation on parameters
    # 
    def __init__(
        self,
        efloat_mode = "float"
    ):
        self.efloat_mode = efloat_mode
        self.efloat = EFloat(efloat_mode)

        self.unique_id_string = ""

        self.section_factory = SectionFactory()
        self.section_list = []


    def setup_pcirexi(
            self,
            contour: Contour,
            interpolation_settings: InterpolationSettings,
            quadrature_settings: QuadratureSettings
    ):
        if self.contour is None or self.interpolation_settings is None or self.contour != contour or \
                self.interpolation_settings != interpolation_settings:
            self.contour = contour
            self.interpolation_settings = interpolation_settings
            self.section_list = self.create_sections()

        self.quadrature_settings = quadrature_settings
        return self._setup_quadrature(quadrature_settings)


    def _setup_quadrature(self, quadrature_settings: QuadratureSettings):
        self.create_coefficients()
        coeffs = REXICoefficients()
        coeffs.alphas = [a for a in self.alphas[:]]
        coeffs.betas = [b for b in self.betas[:]]
        if False: #Wirred hack to get the same oredring as EL-REXI
            coeffs.alphas.reverse()
            coeffs.betas.reverse()
            coeffs.alphas = coeffs.alphas[len(coeffs.alphas)*3//4:]+coeffs.alphas[0:len(coeffs.alphas)*3//4]
            coeffs.betas = coeffs.betas[len(coeffs.alphas)*3//4:]+coeffs.betas[0:len(coeffs.alphas)*3//4]
        coeffs.gamma = 0
        coeffs.unique_id_string = self.getUniqueId()
        coeffs.function_name = quadrature_settings.function_name
        return coeffs


    def getUniqueId(self):
        return "PCIREXI_" + self.unique_id_string


    def create_sections(self):
        if self.interpolation_settings.interp_method.startswith("direct"):
            return self.contour.single_section_list()

        return self.section_factory.generate_sections(
                self.contour.f,
                self.contour.f_derivative,
                self.interpolation_settings.section_number, \
                self.interpolation_settings.interp_per_section,
                self.interpolation_settings.interp_method, \
                equilong=self.interpolation_settings.same_arc_length_for_all_sections,
                sample_points_arc_length_positioned=
                self.interpolation_settings.sample_nodes_based_on_arc_length
            )

    def create_coefficients(self):
        fun = Functions(self.quadrature_settings.function_name, efloat_mode=self.efloat_mode).eval
        cauchy_integrator = CauchyIntegrator(self.section_list, fun, \
                                             self.quadrature_settings.overall_quadrature_points,
                                             self.quadrature_settings.quadrature_method,
                                             arc_length_distribution=self.quadrature_settings.
                                             distribute_quad_points_based_on_arc_length)
        self.alphas = cauchy_integrator.a_s[:]
        self.betas = cauchy_integrator.bs[:]

