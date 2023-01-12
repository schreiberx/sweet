import math

from mule.rexi.pcirexi.contour.Contour import Contour


class RectangleContour(Contour):

    def __init__(self, horizontal_radius: float, vertical_radius: float, center: complex):

        self.center = center
        self.height = vertical_radius
        self.width = horizontal_radius

    # f: [0, 2pi]->complex
    def f(self, alpha: float) -> complex:
        alpha %= (2 * math.pi)
        if alpha <= math.pi / 2:
            return 1j * (-self.height / 2 + alpha * self.height / (math.pi / 2)) + self.width / 2 + self.center
        elif alpha <= math.pi:
            alpha -= math.pi / 2
            return (self.width / 2 - alpha * self.width / (math.pi / 2)) + self.height / 2 * 1j + self.center
        elif alpha <= 3 * math.pi / 2:
            alpha -= math.pi
            return 1j * (self.height / 2 - alpha * self.height / (math.pi / 2)) - self.width / 2 + self.center
        else:
            alpha -= 3 * math.pi / 2
            return (-self.width / 2 + alpha * self.width / (math.pi / 2)) - self.height / 2 * 1j + self.center

    # f_derivative: [0, 2pi]->complex
    def f_derivative(self, alpha: float) -> complex:
        alpha %= (2 * math.pi)
        if alpha <= math.pi / 2:
            return 1j * (self.height / (math.pi / 2))
        elif alpha <= math.pi:
            alpha -= math.pi / 2
            return (-self.width / (math.pi / 2))
        elif alpha <= 3 * math.pi / 2:
            alpha -= math.pi
            return 1j * (-self.height / (math.pi / 2))
        else:
            alpha -= 3 * math.pi / 2
            return (self.width / (math.pi / 2))
