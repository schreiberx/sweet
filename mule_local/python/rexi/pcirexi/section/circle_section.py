from cmath import exp
from math import pi

from mule.rexi.pcirexi.section.section import Section


class CircleSection(Section[float, float]):
    end_a: float

    def __init__(self, start_a, end_a, r=10):
        self.start_a = start_a
        self.end_a = end_a
        self.r = r

    def interpolate(self, a) -> complex:
        return self.r * exp(1j * (self.start_a + (self.end_a - self.start_a) * a))  # +1

    def evaluate_input(self, a) -> complex:
        return 2j * (pi * a / self.end_a)

    def evaluateDerivative(self, a) -> float:
        return self.r * 1j * (self.end_a - self.start_a) * exp(1j * (self.start_a + (self.end_a - self.start_a) * a))

    def evaluate_derivative_parameterized(self, a) -> float:
        return self.r * 1j * exp(1j * (self.start_a + (self.end_a - self.start_a) * a))

    def sub_sections(self, amount: int):
        raise
        return [CircleSection(1)]
