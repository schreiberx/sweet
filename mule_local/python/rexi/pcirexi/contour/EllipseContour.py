import math

from mule.rexi.pcirexi.contour.Contour import Contour


class EllipseContour(Contour):

    def __init__(self, horizontal_radius: float, vertical_radius: float, center: complex):
        self.center = center
        self.vertical_radius = vertical_radius
        self.horizontal_radius = horizontal_radius

    # f: [0, 2pi]->complex
    def f(self, alpha: float) -> complex:
        return self.center + self.horizontal_radius * math.cos(alpha) + 1j * self.vertical_radius * math.sin(alpha)

    # f_derivative: [0, 2pi]->complex
    def f_derivative(self, alpha: float) -> complex:
        return -self.horizontal_radius * math.sin(alpha) + 1j * self.vertical_radius * math.cos(alpha)
