import math
from typing import List

import numpy as np
import scipy.special
from numpy import polyfit

from mule_local.rexi.pcirexi.section.section import Section


class PolynomialSection(Section):
    coefficients = []
    derivation_coefficients = []
    second_derivation_coefficients = []
    x_list = []
    y_list = []

    def __init__(self, start_a, end_a, a_list: List = None, y_list: List = None, coefficients=None, methode='newton'):
        self.start_a = start_a
        self.end_a = end_a
        if a_list is None or y_list is None:
            self.coefficients = coefficients
        else:
            self.a_list = a_list
            self.y_list = y_list
            if methode == 'newton':
                newton_coefficients = self._calculate_newton_coefficients(a_list, y_list)
                self._calculate_coefficients(newton_coefficients)
            elif methode == 'polyfit':
                self.coefficients = [polyfit(a_list, y_list, -1 + len(a_list))[::-1]][0]
            else:
                self._calculate_coefficients_via_linear_equation_system(a_list, y_list)
        self._calculate_derivation_coefficients()

    def _calculate_newton_coefficients(self, a_list, y_list):
        newtonTable = []
        newton_coefficients = []
        newtonTable.append(y_list)
        newton_coefficients.append(y_list[0])
        for k in range(1, len(a_list)):
            newtonTable.append([])
            for i in range(0, len(a_list) - k):
                newtonTable[k].append((newtonTable[k - 1][i + 1] - newtonTable[k - 1][i]) / (a_list[i + k] - a_list[i]))
            newton_coefficients.append(newtonTable[k][0])
        return newton_coefficients

    def interpolate_list(self, a_interpol_list):
        y_interpol_list = []
        for cur_a in a_interpol_list:
            y_interpol_list.append(self.interpolate(cur_a))
        return y_interpol_list

    def interpolate(self, a):
        sum = 0j
        exponent = 1
        for i in range(0, len(self.coefficients)):
            sum += exponent * self.coefficients[i]
            exponent = math.pow(a, i + 1)
        return sum

    def evaluateDerivative(self, a):
        sum = 0j
        exponent = 1
        for i in range(0, len(self.derivation_coefficients)):
            sum += exponent * self.derivation_coefficients[i]
            exponent = math.pow(a, i + 1)
        return sum

    def evaluate_derivative_parameterized(self, a):
        sum = 0j
        exponent = 1
        for i in range(0, len(self.derivation_coefficients)):
            sum += exponent * self.derivation_coefficients[i]
            exponent = math.pow(a, i + 1)
        return sum / (self.end_a - self.start_a)

    def evaluate_second_derivative(self, a):
        sum = 0j
        exponent = 1
        for i in range(0, len(self.second_derivation_coefficients)):
            sum += exponent * self.second_derivation_coefficients[i]
            exponent = math.pow(a, i + 1)
        #    exponent *= a
        return sum

    def evaluate_second_derivative_parameterized(self, a):
        sum = 0j
        exponent = 1
        for i in range(0, len(self.second_derivation_coefficients)):
            sum += exponent * self.second_derivation_coefficients[i]
            exponent = math.pow(a, i + 1)
        #    exponent *= a
        return sum / (self.end_a - self.start_a) / (self.end_a - self.start_a)

    def _calculate_coefficients(self, newton_coefficients):
        self.coefficients = [0j] * len(newton_coefficients)
        for i in range(0, len(newton_coefficients)):
            coefficient_part = [0j] * len(newton_coefficients)
            self._recursive_coefficient(coefficient_part, 0, 1, 0, i - 1)
            coefficient_part = list(map(lambda c_p: c_p * newton_coefficients[i], coefficient_part))
            self.coefficients = [x + y for (x, y) in zip(coefficient_part, self.coefficients)]

    def _recursive_coefficient(self, coefficient_part, degree, value, index, max_index):
        if index > max_index:
            coefficient_part[degree] += value
            return
        self._recursive_coefficient(coefficient_part, degree + 1, value, index + 1, max_index)
        self._recursive_coefficient(coefficient_part, degree, -value * self.a_list[index], index + 1, max_index)
        return

    def _calculate_derivation_coefficients(self):
        self.derivation_coefficients = [0j] * (len(self.coefficients) - 1)
        self.second_derivation_coefficients = [0j] * (len(self.coefficients) - 2)
        for i in range(1, len(self.coefficients)):
            self.derivation_coefficients[i - 1] = self.coefficients[i] * i
        for i in range(2, len(self.coefficients)):
            self.second_derivation_coefficients[i - 2] = self.coefficients[i] * i * (i - 1)

    def sub_sections(self, amount):
        ## First calculate seperation points
        borders = [self.start_a] + [self.end_a] * amount
        for i in range(1, amount):
            last_to_end_length = self.arc_length_start_end(borders[i - 1], self.end_a)
            sub_length = last_to_end_length / (amount - i + 1)
            distance = sub_length - last_to_end_length
            jump_size = (self.end_a - borders[i - 1]) / 2
            steps = 1
            while abs(distance) > sub_length / 1000000000000000:
                print("Steps: " + str(steps))
                steps += 1
                if distance > 0:
                    borders[i] += jump_size
                else:
                    borders[i] -= jump_size
                jump_size /= 2
                distance = sub_length - self.arc_length_start_end(borders[i - 1], borders[i])
        subCoefficient = []
        subSections = []
        ### Paremetrization to alpha \in [0,1]
        for i in range(0, amount):
            subCoefficient.append([0.0] * len(self.coefficients))
            ### coefficients from    a_0*1 + a_1*x                         + a_2*x^2
            ###             to       a_0*1 + a_1*(x*(end-start) + start)   + a_2*(x*(end-start) + start)^2
            for n in range(0, len(self.coefficients)):
                for k in range(0, n + 1):
                    binom = scipy.special.binom(n, k)
                    subCoefficient[i][n - k] += binom * self.coefficients[n] * math.pow(borders[i + 1] - borders[i],
                                                                                        n - k) * math.pow(borders[i], k)
            subSections.append(PolynomialSection(0.0, 1.0, coefficients=subCoefficient[i]))
        return subSections

    def _calculate_coefficients_via_linear_equation_system(self, a_list, y_list):
        A = []
        for a in a_list:
            line = []
            for j in range(len(a_list)):
                line.append(math.pow(a, j))
            A.append(line)
        self.coefficients = np.linalg.solve(A, y_list)
