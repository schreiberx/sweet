from typing import TypeVar, Generic, List

import scipy.integrate as integrate

A_T = TypeVar('A_T')
Y_T = TypeVar('Y_T')


class Section(Generic[A_T, Y_T]):
    start_a: A_T = 0
    end_a: A_T = 0

    def interpolate_list(self, a_interpol_list) -> List[Y_T]:
        y_interpol_list: List[Y_T] = []
        for cur_a in a_interpol_list:
            y_interpol_list.append(self.interpolate(cur_a))
        return y_interpol_list

    def interpolate(self, a) -> Y_T:
        pass

    def evaluateDerivative(self, a) -> Y_T:
        pass

    def arc_length_start_end(self, start: A_T, end: Y_T):
        return integrate.quad(lambda x: abs(self.evaluateDerivative(x)), start, end)[0]

    def sub_sections(self, amount: int):
        pass
