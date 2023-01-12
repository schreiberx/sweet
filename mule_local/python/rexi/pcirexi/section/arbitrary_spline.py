from typing import List, Any
import mpmath

from mule.rexi.pcirexi.gauss_cache import GaussCache
from mule.rexi.pcirexi.section.polynomial_section import PolynomialSection


class ArbitrarySpline:
    a_list = []
    y_list = []
    coefficients: List[List[Any]]

    def __init__(self, degree, start_a, end_a, a_list, y_list):
        if len(a_list) < (degree - 2) * 2 or len(a_list) % (degree - 2) != 0:
            raise
        mpmath.mp.prec = 100
        self.degree = degree
        self.start_a = start_a
        self.end_a = end_a
        self.a_list = a_list
        self.y_list = y_list
        self.g_c = GaussCache()
        self._calculate_coefficients()

    def get_spline_nr(self, a):
        jump_size = len(self.a_list) // 2
        current_pos = jump_size
        jump_size //= 2
        if jump_size == 0:
            jump_size = 1
        while True:
            if a < self.a_list[current_pos]:
                current_pos -= jump_size
            elif current_pos == len(self.a_list) - 1 or a < self.a_list[current_pos + 1]:
                return current_pos // (self.degree - 2)
            else:
                current_pos += jump_size
            jump_size //= 2
            if jump_size == 0:
                jump_size = 1

    def interpolate_list(self, a_interpol_list):
        y_interpol_list = []
        for cur_a in a_interpol_list:
            y_interpol_list.append(self.evaluate(cur_a))
        return y_interpol_list

    def evaluate(self, a):
        spline_nr = self.get_spline_nr(a)
        a = (a - self.a_list[spline_nr * (self.degree - 2)]) / self.h(spline_nr * (self.degree - 2))
        sum = 0j
        exponent = 1
        for i in range(0, len(self.coefficients[spline_nr])):
            sum += exponent * self.coefficients[spline_nr][i]
            exponent = mpmath.mpf.pow(a, i + 1)
            # exponent *= a
        return sum

    def get_polynomial_sections(self):
        poly_sections = []
        for i in range(len(self.coefficients)):
            a_0 = self.a_list[i * (self.degree - 2)]
            if i + 1 == len(self.coefficients):
                a_1 = self.end_a
            else:
                a_1 = self.a_list[(i * (self.degree - 2) + (self.degree - 2)) % len(self.a_list)]
            poly_sections.append(PolynomialSection(a_0, a_1, coefficients=self.coefficients[i]))
        return poly_sections

    def _calculate_coefficients(self):
        basis_polynomials = [self.g_c.spline_basis_polynomials(samples=self.degree - 1, n=i, precision=mpmath.mp.prec)
                             for i in range(0, self.degree - 1 + 2)]
        basis_polynomials_coefficients = [b_p.coefficients for b_p in basis_polynomials]
        second_derivatives = [b_p.second_derivation_coefficients for b_p in basis_polynomials]
        t_0_factors = [s_d[0] for s_d in second_derivatives]
        t_1_factors = [sum(s_d) for s_d in second_derivatives]

        A: List[List[int]] = []
        A_prime = []
        b: List[Any] = []
        step = (self.degree - 2)
        size = len(self.a_list) // step
        for i in range(size):
            row = []
            row_prime = []
            for j in range(size):
                pos = (((j + 1) % size) - i) % size
                h_minus = self.h(i * step - step)
                h = self.h(i * step)
                row.append(t_1_factors[-2] / h_minus if pos == 0 else
                           (t_1_factors[-1] / h_minus - t_0_factors[-2] / h if
                            pos == 1 else -t_0_factors[-1] / h if
                           pos == 2 else 0))
                row_prime.append(pos)
            A.append(row)
            A_prime.append(row_prime)
        for i in range(size):
            h_minus = self.h(i * step - step)
            h = self.h(i * step)
            newValue = 0
            for j in range(step + 1):
                newValue += -t_1_factors[j] * self.y_list[i * step - step + j] / h_minus / h_minus + t_0_factors[j] * \
                            self.y_list[(i * step + j) % len(self.y_list)] / h / h
            b.append(newValue)
        derivations = mpmath.lu_solve(A, b)
        self.coefficients = []
        for i in range(size):
            current_coefficients = [0j] * (self.degree + 1)
            if len(current_coefficients) != len(basis_polynomials_coefficients):
                raise
            for j in range(step + 1):
                for h in range(len(basis_polynomials_coefficients[j])):
                    current_coefficients[h] += basis_polynomials_coefficients[j][h] * self.y_list[
                        (i * step + j) % len(self.y_list)]
            for h in range(len(basis_polynomials_coefficients[-2])):
                current_coefficients[h] += basis_polynomials_coefficients[-2][h] * derivations[i] * self.h(i * step)
            for h in range(len(basis_polynomials_coefficients[-1])):
                current_coefficients[h] += basis_polynomials_coefficients[-1][h] * derivations[
                    (i + 1) % len(derivations)] * self.h(i * step)
            current_coefficients = [complex(c_c) for c_c in current_coefficients]

            self.coefficients.append(current_coefficients)
        return

    def h(self, i):
        if i % (self.degree - 2) != 0:
            raise
        i = i % len(self.a_list)
        if i == len(self.a_list) - (self.degree - 2):
            a_plus = self.end_a
            a = self.a_list[i]
        else:
            a_plus = self.a_list[i + (self.degree - 2)]
            a = self.a_list[i]
        return a_plus - a
