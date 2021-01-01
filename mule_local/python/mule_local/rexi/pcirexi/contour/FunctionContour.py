from cmath import *
import sympy as sym

from mule_local.rexi.pcirexi.contour.Contour import Contour


class FunctionContour(Contour):
    _expression: str
    _expression_derivative: str

    ### expression gets evaluated by system. Parameter is "alpha"
    ### unit_circle: expression = "exp(alpha*1j)"
    def __init__(self, expression: str):
        self._expression = expression

        ### symbolic differentiation:
        sym_exp = sym.parsing.sympy_parser.parse_expr(expression)
        sym_exp_diff = sym.diff(sym_exp)
        self._expression_derivative = str(sym_exp_diff).replace("I", "1j")

    # f: [0, 2pi]->complex
    def f(self, alpha: float) -> complex:
        try:
            return eval(self._expression)
        except:
            print("Problem with evaluating the contour function")
            return 0

    # f_derivative: [0, 2pi]->complex
    def f_derivative(self, alpha: float) -> complex:
        try:
            return eval(self._expression_derivative)
        except:
            print("Problem with evaluating contour derivative")
            return 0
        pass
