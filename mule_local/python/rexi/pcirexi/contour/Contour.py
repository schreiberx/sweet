# interface for contours,
# can be interpolated
# or converted into a single section using ????TODO
import math

from ..section.section import Section


class Contour:
    # f: [0, 2pi]->complex
    def f(self, alpha: float) -> complex:
        pass

    # f_derivative: [0, 2pi]->complex
    def f_derivative(self, alpha: float) -> complex:
        pass

    def single_section_list(self):
        this_contour = self

        class DirectSection(Section):
            def interpolate(self, a) -> complex:
                return this_contour.f(2 * math.pi * a)

            def evaluateDerivative(self, a) -> complex:
                return this_contour.f_derivative(2 * math.pi * a) * 2 * math.pi  ## chain-rule

        return [DirectSection()]
