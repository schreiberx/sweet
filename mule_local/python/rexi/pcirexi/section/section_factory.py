import math
from random import uniform

import numpy as np
from typing import Callable

from scipy import integrate
from scipy.special import roots_chebyt


# binary search for sample points that are distributed based on arc_length
from mule.rexi.pcirexi.gauss_cache import GaussCache
from mule.rexi.pcirexi.section.arbitrary_spline import ArbitrarySpline
from mule.rexi.pcirexi.section.circle_section import CircleSection
from mule.rexi.pcirexi.section.cubic_spline import CubicSpline
from mule.rexi.pcirexi.section.polynomial_section import PolynomialSection


def calculate_arc_length_based_function_points(alphas_function_base: [float], start_section: float, end_section: float,
                                               f: Callable, f_derivative=None):
    DERIVATION_STEP = 10 ** -7
    ACCURACY = 10 ** -9
    if f_derivative is None:
        f_derivative = lambda a: (f(a + DERIVATION_STEP) - f(a - DERIVATION_STEP)) / 2 / DERIVATION_STEP
    arc_length = integrate.quad(lambda x: abs(f_derivative(x)), start_section, end_section)[0]
    accumulated_arc_lengths = [arc_length * (a - start_section) / (end_section - start_section) for a in
                               alphas_function_base]
    current_arc_length = 0
    last_sample_point = start_section
    sample_nodes = []
    for current_goal_arc_length in accumulated_arc_lengths:
        current_step = (end_section - last_sample_point) / 2
        current_sample = last_sample_point + current_step
        current_arc_length = integrate.quad(lambda x: abs(f_derivative(x)), start_section, current_sample)[0]
        # while current_arc_length < current_goal_arc_length * (1 - ACCURACY) or current_arc_length >= current_goal_arc_length * (1 + ACCURACY):
        while abs(current_arc_length - current_goal_arc_length) > ACCURACY and current_step != 0:
            current_step /= 2
            if current_arc_length >= current_goal_arc_length:
                current_sample -= current_step
            else:
                current_sample += current_step
            current_arc_length = integrate.quad(lambda x: abs(f_derivative(x)), start_section, current_sample)[0]
        sample_nodes.append(current_sample)
    return [f(s) for s in sample_nodes]


class SectionFactory:
    g_c: GaussCache

    def __init__(self, g_c=GaussCache()):
        self.g_c = g_c

    def generate_sections(self, f: Callable, f_derivative, section_amount, points_per_section=1,
                          interpolation='equidistant',
                          methode='lgs', circle_radius=10, equilong=False,
                          sample_points_arc_length_positioned=False):
        if interpolation.startswith('cubicsplines'):
            section_type = 'cubicsplines'
        elif interpolation.startswith('spline'):
            section_type = 'spline'
        elif interpolation.startswith('circle'):
            section_type = 'circle'
        else:
            section_type = 'polynomials'
            interpolation = interpolation.split(' ')[0]
            interpolation = interpolation.split(',')[0]
        if equilong:
            section_start_list = self.generate_equilong_section_start_list(f, f_derivative, section_amount)
        else:
            section_start_list = np.linspace(0, 2 * np.pi * (1 - 1 / section_amount), section_amount)

        if section_type == 'polynomials':
            interpolation_alphas = [a / 2 + 0.5 for a in
                                    self._alphas_for_interpolation(interpolation, points_per_section)]
            section_list = []
            for n in range(0, section_amount):
                if n == section_amount - 1:
                    end_section = 2 * np.pi
                else:
                    end_section = section_start_list[n + 1]
                alphas_function_base = [section_start_list[n] + (end_section - section_start_list[n]) * a for a in
                                        interpolation_alphas]
                if sample_points_arc_length_positioned:
                    function_points = calculate_arc_length_based_function_points(alphas_function_base,
                                                                                 section_start_list[n], end_section, f, f_derivative)
                else:
                    function_points = [f(a) for a in alphas_function_base]
                section = PolynomialSection(section_start_list[n],
                                            end_section, interpolation_alphas,
                                            function_points, methode=methode)
                section_list.append(section)
            return section_list
        elif section_type == 'cubicsplines':
            interpolation_alphas = np.linspace(0, 2 * math.pi * (1 - 1 / section_amount), section_amount)
            #            interpolation_alphas = [0] + sorted([uniform(0, 2*math.pi*(1-1/section_amount)) for i in range(1, section_amount)])
            sample_alphas_orig = interpolation_alphas
            sections = CubicSpline(0, 2 * math.pi, interpolation_alphas,
                                   [f(a) for a in sample_alphas_orig]).get_polynomial_sections()
            sample_alphas_base = np.linspace(0, 1, 40)
            interpolation_alphas = []
            ys = []
            for section in sections:
                interpolation_alphas += [a * (section.end_a - section.start_a) + section.start_a for a in
                                         sample_alphas_base]
                ys += [section.interpolate(a) for a in sample_alphas_base]
            # gui.plot_list([y.real for y in ys], [y.imag for y in ys])
            return sections
        elif section_type == 'spline':
            interpolation_alphas = np.linspace(0, 2 * math.pi * (1 - 1 / section_amount / (points_per_section - 1)),
                                               section_amount * (points_per_section - 1))
            #            interpolation_alphas = [0] + sorted([uniform(0, 2*math.pi*(1-1/section_amount)) for i in range(1, section_amount)])
            sample_alphas_orig = interpolation_alphas
            spline = ArbitrarySpline(points_per_section + 1, 0, 2 * math.pi, interpolation_alphas,
                                     [f(a) for a in sample_alphas_orig])
            sections = spline.get_polynomial_sections()
            if methode == 'basis_condition_return':
                return spline.basis_condition
            if methode == 'overall_condition_return':
                return spline.overall_condition
            sample_alphas_base = np.linspace(0, 1, 40)
            interpolation_alphas = []
            ys = []
            for section in sections:
                interpolation_alphas += [a * (section.end_a - section.start_a) + section.start_a for a in
                                         sample_alphas_base]
                ys += [section.interpolate(a) for a in sample_alphas_base]
            # gui.plot_list([y.real for y in ys], [y.imag for y in ys])
            return sections
        elif section_type == 'circle':
            section_list = []
            interpolation_alphas = np.linspace(0, 2 * math.pi * (1), section_amount + 1)
            # interpolation_alphas = [0] + sorted([uniform(0, 2*math.pi*(1-1/section_amount)) for i in range(1, section_amount)]) + [2*math.pi]
            for n in range(0, section_amount):
                section_list.append(CircleSection(interpolation_alphas[n], interpolation_alphas[n + 1], circle_radius))
            sample_alphas_base = np.linspace(0, 1, 40)
            interpolation_alphas = []
            ys = []
            for section in section_list:
                interpolation_alphas += [a * (section.end_a - section.start_a) + section.start_a for a in
                                         sample_alphas_base]
                ys += [section.evaluate_derivative_parameterized(a) for a in sample_alphas_base]
            # gui.plot_list(interpolation_alphas, ys)
            return section_list

    def generate_equilong_section_start_list(self, f, f_derivative, section_amount):
        section_start_list = [0] * section_amount
        accuracy = 10 ** -10
        left_length = \
            integrate.quad(lambda x: abs(f_derivative(x)), 0, 2 * np.pi)[0]
        current_pos = 0
        current_width = 2 * np.pi
        step = current_width / 2
        for i in range(1, section_amount):
            goal_length = left_length / (section_amount - i + 1)
            cur_length = \
                integrate.quad(lambda x: abs(f_derivative(x)), current_pos, current_pos + current_width)[0]
            while cur_length < goal_length * (1 - accuracy) or cur_length >= goal_length * (1 + accuracy):
                if cur_length >= goal_length * (1 + accuracy):
                    current_width -= step
                    step /= 2
                else:
                    current_width += step
                    step /= 2
                cur_length = \
                    integrate.quad(lambda x: abs(f_derivative(x)),
                                   current_pos, current_pos + current_width)[0]
            current_pos = current_pos + current_width
            section_start_list[i] = current_pos
            left_length -= cur_length
            current_width = left_length
            step = current_width / 2
        return section_start_list

    def _alphas_for_interpolation(self, interpolation, sample_points):
        sample_alphas = []
        if interpolation == 'equidistant':
            # print("Ip ed")
            # with 1/(2*section) distance to bodrders
            # sample_alphas = np.linspace(-1 + 1/(samples_per_section*2), 1 - 1/(2*samples_per_section), samples_per_section)
            sample_alphas = np.linspace(-1, 1, sample_points)
        elif interpolation == 'lobatto':
            # print("Ip lobatto")
            base_nodes, _ = self.g_c.gauss_lobatto(sample_points, 20)
            sample_alphas = [float(b_n) for b_n in base_nodes]
        elif interpolation == 'legendre':
            # print("Ip legendre")
            base_nodes, _ = self.g_c.gauss_legendre(sample_points, 20)
            sample_alphas = [float(b_n) for b_n in base_nodes]
        elif interpolation == 'chebychev':
            # print("Ip chebychev")
            base_nodes, _ = roots_chebyt(sample_points)
            sample_alphas = base_nodes
        return sample_alphas
