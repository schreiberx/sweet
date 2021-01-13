#! /usr/bin/env python3
#
# Author: Raphael Schilling <raphael.schilling@tum.de>
#

from mule_local.rexi.pcirexi.PCIREXI import PCIREXI
from mule_local.rexi.pcirexi.contour.RectangleContour import RectangleContour
from mule_local.rexi.pcirexi.section.InterpolationSettings import InterpolationSettings
from mule_local.rexi.pcirexi.section.QuadratureSettings import QuadratureSettings
from mule_local.rexi.EFloat import *



class LRREXI:

    def __init__(
        self,
        efloat_mode = "float"
    ):
        self.efloat_mode = efloat_mode
        self.efloat = EFloat(efloat_mode)

        self.unique_id_string = ""


    def setup(
            self,
            function_name:str,
            width:float,
            height:float,
            center:complex,
            N:int
    ):
        pcirexi = PCIREXI()
        contour = RectangleContour(width, height, center)
        i_s = InterpolationSettings(4, 2, 'equidistant', False, False)
        q_s = QuadratureSettings(overall_quadrature_points=N,
                                 quadrature_method="legendre",
                                 distribute_quad_points_based_on_arc_length=True,
                                 function_name=function_name)
        coeffs = pcirexi.setup_pcirexi(contour=contour,
                                            interpolation_settings=i_s,
                                            quadrature_settings=q_s)

        self.unique_id_string = function_name
        self.unique_id_string += "_N" + str(N)
        self.unique_id_string += "_w" + str(width)
        self.unique_id_string += "_h" + str(height)
        self.unique_id_string += "_c" + str(center)

        coeffs.unique_id_string = self.unique_id_string
        return coeffs


    def getUniqueId(self):
    	return "LRREXI_"+self.unique_id_string

