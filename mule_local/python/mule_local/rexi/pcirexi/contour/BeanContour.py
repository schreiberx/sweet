import math

from mule_local.rexi.pcirexi.contour.Contour import Contour


class BeanContour(Contour):

    def __init__(self, horizontal_radius: float, vertical_radius: float, center: float, compression=0.25):
        self.compression = compression
        self.center = center
        self.vertical_radius = vertical_radius
        self.horizontal_radius = horizontal_radius

    # f: [0, 2pi]->complex
    def f(self, alpha: float) -> complex:
        return 1j * math.sin(alpha) * self.vertical_radius + math.cos(alpha) * (
                1 - math.pow(math.cos(alpha), 3) * self.compression) * self.horizontal_radius + self.center

    # f_derivative: [0, 2pi]->complex
    def f_derivative(self, alpha: float) -> complex:
        return 4 * self.horizontal_radius * self.compression * math.sin(alpha) * math.pow(math.cos(alpha),
                                                                                          3) + 1j * self.vertical_radius * math.cos(
            alpha) - self.horizontal_radius * math.sin(alpha)
