import math
from typing import List, Any
import numpy as np

from mule_local.rexi.pcirexi.section.polynomial_section import PolynomialSection


class CubicSpline:
    a_list = []
    y_list = []
    coefficients: List[List[Any]]

    def __init__(self, start_a, end_a, a_list, y_list):
        self.start_a = start_a
        self.end_a = end_a
        self.a_list = a_list
        self.y_list = y_list
        self.calculate_coefficients()
        self.calculate_derivation_coefficients()

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
                return current_pos
            else:
                current_pos += jump_size
            jump_size //= 2
            if jump_size == 0:
                jump_size = 1

    def interpolate(self, a_interpol_list):
        y_interpol_list = []
        for cur_a in a_interpol_list:
            y_interpol_list.append(self.evaluatePolynome(cur_a))
        return y_interpol_list

    def evaluate(self, a):
        spline_nr = self.get_spline_nr(a)
        a = (a - self.a_list[spline_nr]) / self.h(spline_nr)
        sum = 0j
        exponent = 1
        for i in range(0, len(self.coefficients[spline_nr])):
            sum += exponent * self.coefficients[spline_nr][i]
            exponent = math.pow(a, i + 1)
            # exponent *= a
        return sum

    def evaluateDerivative(self, a):
        spline_nr = self.get_spline_nr(a)
        a = (a - self.a_list[spline_nr]) / self.h(spline_nr)
        sum = 0j
        exponent = 1
        for i in range(0, len(self.derivation_coefficients[spline_nr])):
            sum += exponent * self.derivation_coefficients[spline_nr][i]
            exponent = math.pow(a, i + 1)
        return sum / self.h(spline_nr)

    def calculateDerivationCoefficients(self):
        self.derivation_coefficients = [0j] * (len(self.coefficients) - 1)
        for i in range(1, len(self.coefficients)):
            self.derivation_coefficients[i - 1] = self.coefficients[i] * i

    def get_polynomial_sections(self):
        poly_sections = []
        for i in range(len(self.a_list)):
            a_0 = self.a_list[i]
            if i + 1 == len(self.a_list):
                a_1 = self.end_a
            else:
                a_1 = self.a_list[(i + 1) % len(self.a_list)]
            poly_sections.append(PolynomialSection(a_0, a_1, coefficients=self.coefficients[i]))
        return poly_sections

    def calculate_coefficients(self):
        A: List[List[int]] = []
        b: List[Any] = []
        for i in range(len(self.a_list)):
            row = []
            for j in range(len(self.a_list)):
                pos = (((j + 1) % len(self.a_list)) - i) % len(self.a_list)
                h_minus = self.h(i - 1)
                h = self.h(i)
                row.append(1 / h_minus if pos == 0 else
                           (2 / h_minus + 2 / h if
                            pos == 1 else 1 / h if
                           pos == 2 else 0))
            A.append(row)
        for i in range(len(self.a_list)):
            h_minus = self.h(i - 1)
            h = self.h(i)
            b.append(3 * ((self.y_list[i] - self.y_list[(i - 1) % len(self.a_list)]) / (h_minus * h_minus)
                          + (self.y_list[(i + 1) % len(self.a_list)] - self.y_list[i]) / (h * h)))
        #            b.append((self.y_list[(i + 1) % len(self.a_list)] - self.y_list[(i-1) % len(self.a_list)])/(self.a_list[1]-self.a_list[0])*3)
        derivations = np.linalg.solve(A, b)
        self.coefficients = []
        for i in range(len(self.a_list)):
            current_coefficients = [0j] * 4
            current_coefficients[0] = self.y_list[i]
            h = self.h(i)
            current_coefficients[1] = derivations[i] * h
            current_coefficients[2] = -3 * self.y_list[i] + 3 * self.y_list[(i + 1) % len(self.a_list)] - 2 * \
                                      derivations[i] * h - derivations[(i + 1) % len(self.a_list)] * h
            current_coefficients[3] = 2 * self.y_list[i] - 2 * self.y_list[(i + 1) % len(self.a_list)] + derivations[
                i] * h + derivations[(i + 1) % len(self.a_list)] * h
            self.coefficients.append(current_coefficients)

    def calculate_coefficients2(self):
        A: List[List[int]] = []
        b: List[Any] = []
        for i in range(len(self.a_list)):
            row = []
            for j in range(len(self.a_list)):
                pos = (((j + 1) % len(self.a_list)) - i) % len(self.a_list)
                row.append(1 / 6 if pos == 0 else
                           (2 / 3 if
                            pos == 1 else 1 / 6 if
                           pos == 2 else 0))
            A.append(row)
            print(row)

        for i in range(len(self.a_list)):
            h_minus = self.h(i - 1)
            h = self.h(i)
            y_i_p = self.y_list[(i + 1) % len(self.a_list)]
            y_i = self.y_list[i]
            y_i_m = self.y_list[(i - 1) % len(self.a_list)]
            b.append(y_i_p - 2 * y_i + y_i_m)
        #            b.append((self.y_list[(i + 1) % len(self.a_list)] - self.y_list[(i-1) % len(self.a_list)])/(self.a_list[1]-self.a_list[0])*3)
        derivations = np.linalg.solve(A, b)
        self.coefficients = []
        for i in range(len(self.a_list)):
            current_coefficients = [0j] * 4
            current_coefficients[0] = self.y_list[i]
            y_i_p = self.y_list[(i + 1) % len(self.y_list)]
            y_i = self.y_list[i]
            M_i = derivations[i]
            M_i_p = derivations[(i + 1) % len(derivations)]
            current_coefficients[1] = -1 / 2 * M_i + y_i_p - y_i - (M_i_p - M_i) / 6
            current_coefficients[2] = 1 / 2 * M_i
            current_coefficients[3] = -1 / 6 * M_i + 1 / 6 * M_i_p
            self.coefficients.append(current_coefficients)

    def calculate_coefficients3(self):
        A: List[List[int]] = []
        b: List[Any] = []
        for i in range(len(self.a_list)):
            row = []
            for j in range(len(self.a_list)):
                pos = (((j + 1) % len(self.a_list)) - i) % len(self.a_list)
                h_minus = self.h(i - 1)
                h = self.h(i)
                row.append(h_minus / 6 if pos == 0 else
                           ((h_minus + h) / 3 if
                            pos == 1 else h / 6 if
                           pos == 2 else 0))
            A.append(row)
            print(row)

        for i in range(len(self.a_list)):
            h_minus = self.h(i - 1)
            h = self.h(i)
            y_i_p = self.y_list[(i + 1) % len(self.a_list)]
            y_i = self.y_list[i]
            y_i_m = self.y_list[(i - 1) % len(self.a_list)]
            b.append((y_i_p - y_i) / h - (y_i - y_i_m / h_minus))
        #            b.append((self.y_list[(i + 1) % len(self.a_list)] - self.y_list[(i-1) % len(self.a_list)])/(self.a_list[1]-self.a_list[0])*3)
        derivations = np.linalg.solve(A, b)
        self.coefficients = []
        for i in range(len(self.a_list)):
            current_coefficients = [0j] * 4
            h = self.h(i)
            y_i_p = self.y_list[(i + 1) % len(self.y_list)]
            y_i = self.y_list[i]
            M_i = derivations[i]
            M_i_p = derivations[(i + 1) % len(derivations)]
            # current_coefficients[0] = self.y_list[i]
            # current_coefficients[1] = -1 / 2 *h**2 *M_i + y_i_p - y_i - h**2*(M_i_p - M_i)/6
            # current_coefficients[2] = 1/2*M_i*h**2
            # current_coefficients[3] = -1/6*M_i*h**2 + 1/6*M_i_p*h**2
            current_coefficients[0] = self.y_list[i]
            current_coefficients[1] = -1 / 2 * h ** 2 * M_i + y_i_p - y_i - (h ** 2) * (M_i_p - M_i) / 6
            current_coefficients[2] = 1 / 2 * h ** 2 * M_i
            current_coefficients[3] = -1 / 6 * h ** 2 * M_i + 1 / 6 * h ** 2 * M_i_p
            self.coefficients.append(current_coefficients)

    def calculate_coefficients4(
            self):  # setting s_n'(alpha_n + 0*(alpha_(n+1)-alpha_n)) =s_n_1'(alpha_n + 1*(alpha_(n+1)-alpha_n))
        A: List[List[int]] = []
        b: List[Any] = []
        for i in range(len(self.a_list)):
            row = []
            for j in range(len(self.a_list)):
                pos = (((j + 1) % len(self.a_list)) - i) % len(self.a_list)
                h_minus = self.h(i - 1)
                h = self.h(i)
                row.append((h_minus ** 2) / 6 if pos == 0 else
                           ((h_minus ** 2) / 3 + (h ** 2) / 3 if
                            pos == 1 else (h ** 2) / 6 if
                           pos == 2 else 0))
            A.append(row)
            print(row)

        for i in range(len(self.a_list)):
            h_minus = self.h(i - 1)
            h = self.h(i)
            y_i_p = self.y_list[(i + 1) % len(self.a_list)]
            y_i = self.y_list[i]
            y_i_m = self.y_list[(i - 1) % len(self.a_list)]
            b.append((y_i_p - y_i) - (y_i - y_i_m))
        #            b.append((self.y_list[(i + 1) % len(self.a_list)] - self.y_list[(i-1) % len(self.a_list)])/(self.a_list[1]-self.a_list[0])*3)
        derivations = np.linalg.solve(A, b)
        self.coefficients = []
        for i in range(len(self.a_list)):
            current_coefficients = [0j] * 4
            current_coefficients[0] = self.y_list[i]
            h = self.h(i)
            y_i_p = self.y_list[(i + 1) % len(self.y_list)]
            y_i = self.y_list[i]
            M_i = derivations[i]
            M_i_p = derivations[(i + 1) % len(derivations)]
            current_coefficients[1] = -1 / 2 * h ** 2 * M_i + y_i_p - y_i - h ** 2 * (M_i_p - M_i) / 6
            current_coefficients[2] = 1 / 2 * M_i * h ** 2
            current_coefficients[3] = -1 / 6 * M_i * h ** 2 + 1 / 6 * M_i_p * h ** 2
            self.coefficients.append(current_coefficients)

    def h(self, i):
        i = i % len(self.a_list)
        if i == len(self.a_list) - 1:
            a_plus = self.end_a
            a = self.a_list[i]
        else:
            a_plus = self.a_list[i + 1]
            a = self.a_list[i]
        return a_plus - a

    def calculate_derivation_coefficients(self):
        self.derivation_coefficients = []
        for i in range(len(self.coefficients)):
            current_derivation_coefficients = [0j] * (len(self.coefficients[0]) - 1)
            for j in range(1, len(self.coefficients[0])):
                current_derivation_coefficients[j - 1] = self.coefficients[i][j] * j
            self.derivation_coefficients.append(current_derivation_coefficients)
