import cmath
import numpy as np


def ellipse(radius, rfactor, ifactor):
    def contour(alpha, to_string=False):
        if to_string:
            return "ellipse(" + str(radius * rfactor) + ", " + str(radius * ifactor) + ")"
        value = radius * cmath.exp(1j * alpha)
        return value.real * rfactor + value.imag * ifactor * 1j

    return contour


def stick(radius, cube_radius):
    def contour(alpha, to_string=False):
        if to_string:
            return "stick(" + str(radius) + ", " + str(cube_radius) + ")"
        value = radius * cmath.exp(1j * alpha)
        if alpha < np.pi / 4:
            return radius + 1j * alpha / (np.pi / 4) * cube_radius
        elif alpha < 3 / 4 * np.pi:
            return radius * cmath.exp((alpha - np.pi / 4) * 2j) + cube_radius * 1j
        elif alpha < 5 / 4 * np.pi:
            return 0 - radius + 1j * (np.pi - alpha) / (np.pi / 4) * cube_radius
        elif alpha < 7 / 4 * np.pi:
            return radius * cmath.exp(np.pi * 1j + (alpha - np.pi * 5 / 4) * 2j) - cube_radius * 1j
        else:
            return radius + 1j * (alpha - 2 * np.pi) / (np.pi / 4) * cube_radius

    return contour


def circle(radius):
    def contour(alpha, to_string=False):
        if to_string:
            return "circle(" + str(radius) + ")"
        return 1j * radius * cmath.exp(alpha * 1j)

    return contour
