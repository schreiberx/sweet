from functools import reduce, lru_cache
import numpy as np
from typing import List, Any, Callable
import scipy.integrate as integrate

from mule_local.rexi.pcirexi.gauss_cache import GaussCache
from mule_local.rexi.pcirexi.section import section


def _complex_quad(func, a, b):
    real = integrate.quad((lambda a: func(a).real), a, b)[0]
    imag = integrate.quad((lambda a: func(a).imag), a, b)[0]
    return real + 1j * imag

class CauchyIntegrator:
    target_function: Callable[[Any], Any]  # exp for REXI
    bs: List[complex]  # numerator of terms
    a_s: List[complex]  # part of the denominator
    terms_number: int  # overall number of REXI terms
    terms_section_list: int  # REXI terms per section
    sections: List[section.Section]  # list of contour functions

    def __init__(self, sections, target_function, terms_number=20, integral='trapeze', arc_length_distribution=False,
                 normalize_for_zero=False, g_c=GaussCache()):
        self.g_c = g_c
        self.target_function = target_function
        self.sections = sections
        self.terms_number = terms_number
        self.integral = integral
        self._calculate_quad_points_for_sections(arc_length_distribution)
        if integral.startswith("trapeze") or integral.startswith("lobatto"):
            self.edge_to_edge = True
            # Temporarly increase terms. Later undone by merging values from bs and a_s
            self.terms_section_list = [t + 1 for t in self.terms_section_list]
        else:
            self.edge_to_edge = False
        self._calculateTable()
        if normalize_for_zero:
            self._normalize_for_zero()

    @lru_cache(maxsize=32)
    def get_quadrature_points_and_weights(self, g_c, integral, terms_number):
        # if terms_number % self.terms_per_section != 0:
        #    print("#terms is not dividable by operations_per_section")
        if (not integral.startswith('trapeze')) and (not integral.startswith('midpoint')):
            if integral.startswith('lobatto'):
                if terms_number < 2:
                    print("Lobatto needs at least 2 base_nodes")
                    raise
                print("Ig: lobatto")
                base_nodes, weights_original = g_c.gauss_lobatto(terms_number, 20)
            elif integral.startswith('legendre'):
                base_nodes, weights_original = g_c.gauss_legendre(terms_number, 20)
            elif integral.startswith('chebychev'):
                raise
            else:
                print("Unknown Interation")
                raise
            weights = [float(w) / 2 for w in weights_original]
            base_nodes = [float(b) / 2 + 0.5 for b in base_nodes]
        elif integral.startswith('trapeze'):
            if terms_number == 0:
                base_nodes = []
                weights = []
            else:
                base_nodes = np.linspace(0, 1, terms_number)
                if terms_number <= 2:
                    weights = [1 / terms_number] * terms_number
                else:
                    devisor = (2 * (terms_number - 1))
                    weights = [1 / devisor] + [2 / devisor] * (terms_number - 2) + [1 / devisor]
        elif integral.startswith('midpoint'):
            if terms_number == 0:
                base_nodes = []
                weights = []
            else:
                base_nodes = np.linspace(1 / (2 * terms_number), 1 - 1 / (2 * terms_number),
                                         terms_number)
                weights = [1 / (terms_number)] * terms_number

        else:
            raise
        return (base_nodes, weights)

    def _calculateTable(self):
        # calculates numerator (bs) and denominator addend (as)
        self.bs = []
        self.a_s = []
        for j in range(len(self.sections)):
            current_section = self.sections[j]
            terms = self.terms_section_list[j]
            base_nodes, weights = self.get_quadrature_points_and_weights(self.g_c, self.integral, terms)
            for i in range(terms):
                # prepare one REXI term
                alpha_for_current_term = base_nodes[i]
                contour_pos = current_section.interpolate(alpha_for_current_term)
                contour_derivative = current_section.evaluateDerivative(alpha_for_current_term)
                function_evaluation_at_contour_pos = self.target_function(contour_pos)
                b = -1 / (2j * np.pi) * function_evaluation_at_contour_pos * contour_derivative * weights[i]
                self.bs.append(b)
                self.a_s.append(-contour_pos)
        if self.edge_to_edge:
            # Undo temporary increase of terms
            current_transition = 0
            for i in range(len(self.sections)):
                # Merge values at equal contour position
                self.bs[current_transition] += self.bs[current_transition - 1]
                self.a_s[current_transition] = self.a_s[current_transition] / 2 + self.a_s[current_transition - 1] / 2
                current_transition += self.terms_section_list[i]
            current_unwanted = 0
            for i in range(len(self.sections)):
                # Pop unwanted values
                current_unwanted += self.terms_section_list[i] - 1
                self.bs.pop(current_unwanted)
                self.a_s.pop(current_unwanted)

    def _normalize_for_zero(self):
        current = self.approximate_target_function(0)
        actual = self.target_function(0)
        factor = actual.real / current.real
        self.bs = [b * factor for b in self.bs]

    def approximate_target_function(self, x):
        sum: complex = 0j
        for s in range(0, len(self.bs)):
            sum += self.bs[s] / (self.a_s[s] + x)
        return sum

    def approximate_target_function_using_scipy_quad(self, x):
        sections_sum = 0j
        for current_section in self.sections:
            def cauchy_integrand(alpha):
                contour_pos = current_section.interpolate(alpha)
                contour_derivative = current_section.evaluateDerivative(alpha)
                return self.target_function(contour_pos) * contour_derivative / (contour_pos - x)

            sections_sum += _complex_quad(cauchy_integrand, 0, 1)
        return sections_sum / (2j * np.pi)

    def _get_section(self, a):
        jump_size = len(self.sections) // 2
        current_pos = jump_size
        jump_size //= 2
        if jump_size == 0:
            jump_size = 1
        while True:
            if a < self.sections[current_pos].start_a:
                current_pos -= jump_size
            elif current_pos == len(self.sections) - 1 or a < self.sections[current_pos].end_a:
                return current_pos
            else:
                current_pos += jump_size
            jump_size //= 2
            if jump_size == 0:
                jump_size = 1

    def calc_max_error_in_interval(self, lower_i=-8, higher_i=8, samples=1000):
        values = np.linspace(lower_i * 1j, higher_i * 1j, samples)
        deviations = [abs(self.target_function(a) - self.approximate_target_function(a)) for a in values]
        max_deviation = max(deviations)
        return max_deviation

    def calc_max_error_in_intervall_via_scipy_quad(self, lower_i=-8, higher_i=8, samples=1000):
        values = np.linspace(lower_i * 1j, higher_i * 1j, samples)
        deviations = [abs(self.target_function(a) - self.approximate_target_function_using_scipy_quad(a)) for a in
                      values]
        max_deviation = max(deviations)
        return max_deviation

    def _calculate_quad_points_for_sections(self, arc_length_distribution):
        if not arc_length_distribution:
            self.terms_section_list = [self.terms_number // len(self.sections)] * len(self.sections)
            return
        self.terms_section_list = [0] * len(self.sections)
        section_lengths = [s.arc_length_start_end(0, 1) for s in self.sections]
        contour_length = reduce(float.__add__, section_lengths, 0.0)
        length_per_quad_point = contour_length / self.terms_number
        self.terms_section_list = [int(l / length_per_quad_point) for l in section_lengths]
        current_terms = reduce(int.__add__, self.terms_section_list, 0)
        for _ in range(0, self.terms_number - current_terms):
            max_value = -1
            max_i_so_far = 0
            for i in range(0, len(section_lengths)):
                current_length = section_lengths[i] - self.terms_section_list[i] * length_per_quad_point
                if max_value < current_length:
                    max_i_so_far = i
                    max_value = current_length
            self.terms_section_list[max_i_so_far] += 1
