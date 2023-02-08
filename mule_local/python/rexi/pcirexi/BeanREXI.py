#! /usr/bin/env python3
#
# Author: Raphael Schilling <raphael.schilling@tum.de>
#

from mule.rexi.pcirexi.PCIREXI import PCIREXI
from mule.rexi.pcirexi.contour.BeanContour import BeanContour
from mule.rexi.pcirexi.contour.RectangleContour import RectangleContour
from mule.rexi.pcirexi.section.InterpolationSettings import InterpolationSettings
from mule.rexi.pcirexi.section.QuadratureSettings import QuadratureSettings
from mule.rexi.EFloat import *


class BeanREXI:
    def __init__(
        self,
        efloat_mode = "float"
    ):
        self.efloat_mode = efloat_mode
        self.efloat = EFloat(efloat_mode)

        self.unique_id_string = ""


    def setup(self, function_name, horizontal_radius:float, vertical_radius:float, center:complex, N):
        pcirexi = PCIREXI()
        contour = BeanContour(horizontal_radius, vertical_radius, center, compression=0.25)
        i_s = InterpolationSettings(1, 0, 'direct', False, False)
        q_s = QuadratureSettings(overall_quadrature_points=N,
                                 quadrature_method="midpoint",
                                 distribute_quad_points_based_on_arc_length=False,
                                 function_name=function_name)
        coeffs = pcirexi.setup_pcirexi(contour=contour,
                                            interpolation_settings=i_s,
                                            quadrature_settings=q_s)

        unique_id_string = "BeanREXI_"
        unique_id_string += function_name
        unique_id_string += "_N" + str(N)
        unique_id_string += "_hr" + str(horizontal_radius)
        unique_id_string += "_vr" + str(vertical_radius)
        unique_id_string += "_c" + str(center)

        coeffs.unique_id_string = unique_id_string
        return coeffs
