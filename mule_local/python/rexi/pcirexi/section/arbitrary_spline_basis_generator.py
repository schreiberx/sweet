import mpmath

from mule.rexi.pcirexi.section.polynomial_section import PolynomialSection


def generate_basis_polynomials(samples, n, precision, calc_condition=False):
    mpmath.mp.prec = precision
    A = []
    for i in range(samples):
        line = []
        for j in range(samples + 2):
            t = mpmath.mpf(1) / (samples - 1) * i
            line.append(t ** j)
        A.append(line)
    for t in range(2):
        line = []
        for j in range(samples + 2):
            if j == 0:
                line.append(mpmath.mpf(0.0))
            else:
                line.append(mpmath.mpf(j) * (t ** (j - 1)))
        A.append(line)
    b = [mpmath.mpf(0)] * n + [mpmath.mpf(1)] + [mpmath.mpf(0)] * (samples + 2 - n - 1)
    if calc_condition:
        basis_condition = mpmath.cond(A)
        return basis_condition
    x = mpmath.lu_solve(A, b)
    # x = [float(v) for v in mpmath.mp.lu_solve(A,b)]
    return PolynomialSection(0, 1, coefficients=[element for element in x])
