#! /usr/bin/env python3

import numpy as np
import sys


def _quad_coeffs_hq(M, collocation_type, digits=20):
    if M == 1:
        x = np.array([0.0])
        w = np.array([2.0])

    elif collocation_type == "gauss_legendre":
        from sympy.integrals.quadrature import gauss_legendre
        x, w = gauss_legendre(M, digits)

    elif collocation_type == "gauss_lobatto":
        from sympy.integrals.quadrature import gauss_lobatto
        x, w = gauss_lobatto(M, digits)

    elif collocation_type == "gauss_hermite":
        from sympy.integrals.quadrature import gauss_hermite
        x, w = gauss_hermite(M, 30)
        x = np.array(x, dtype=float)
        w = np.array(w, dtype=float)


    elif collocation_type == "gauss_jacobi":
        from sympy.integrals.quadrature import gauss_jacobi
        x, w = gauss_jacobi(M, 0, 0, 30)
        x = np.array(x, dtype=float)
        w = np.array(w, dtype=float)


    elif collocation_type == "gauss_chebyshev_u":
        from sympy.integrals.quadrature import gauss_chebyshev_u
        x, w = gauss_chebyshev_u(M, 30)
        x = np.array(x, dtype=float)
        w = np.array(w, dtype=float)


    elif collocation_type == "gauss_chebyshev_t":
        from sympy.integrals.quadrature import gauss_chebyshev_t
        x, w = gauss_chebyshev_t(M, 30)
        x = np.array(x, dtype=float)
        w = np.array(w, dtype=float)


    else:
        raise Exception("Unknown collocation method '"+str(collocation_type)+"'")

    assert len(x) == M

    return x, w



def quad_coeffs(M, collocation_type):

    if M == 1:
        x = np.array([0.0])
        w = np.array([2.0])

    elif collocation_type == "equidistant":
        x = np.linspace(0, 1, M, endpoint=True)
        w = np.ones(M)/(M-1)
        w[0] *= 0.5
        w[-1] *= 0.5

    elif collocation_type == "geometric":
        x = np.zeros(M)
        for i in range(M//2):
            x[i] = -(2.0**-(i+1))
            x[M-1-i] = (2.0**-(i+1))
        w = None

    elif collocation_type in ["gauss_legendre", "gauss_lobatto", "gauss_hermite", "gauss_jacobi", "gauss_chebyshev_u", "gauss_chebyshev_t"]:
        x, w = _quad_coeffs_hq(M, collocation_type, digits=20)
        x = np.array(x, dtype=float)
        w = np.array(w, dtype=float)

    elif collocation_type == "chebyshev_gauss_lobatto":
        p = np.linspace(0, 1, M, endpoint=True)
        x = -np.cos(np.pi*p)
        w = np.zeros(M, dtype=float)
        w[:] = 1.0/(M+1)
        w[0] *= 0.5
        w[-1] *= 0.5

    else:
        raise Exception("Unknown collocation method '"+str(collocation_type)+"'")

    assert len(x) == M

    return x, w

