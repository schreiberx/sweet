#! /usr/bin/env python3

"""
Copyright 2017, Jed Brown

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
Author: Jed Brown

Other contributions by: Matthew Normile, Martin Schreiber
"""


import numpy
import mule_local.rexi.EFloat as ef


floatmode = "float"
#floatmode = "mpfloat"



def _golubwelsch(n, lo=-1, hi=1, radau=False):

    if floatmode == "float":

        beta = .5 / numpy.sqrt(1-(2*(numpy.arange(1,n)))**(-2.)) # 3-term recurrence coeffs
        T = numpy.diag(beta,1) + numpy.diag(beta,-1);     # Jacobi matrix
        if radau:
            # Eq. 2.4 of Gautschi, Gauss–Radau formulae for Jacobi and Laguerre weight functions, 2000
            if n == 1:
                T[-1,-1] = 1
            else:
                T[-1,-1] = 1 - 2*(n-1)**2 / (2*(n-1)*(2*(n-1)+1))

        D, V = numpy.linalg.eigh(T)                      # Eigenvalue decomposition

        i = numpy.argsort(D)         # Legendre points
        x = D[i]
        w = 2*V[0,i]**2                   # Quadrature weights
        x = 0.5*((hi+lo) + (hi-lo)*x)
        w = 0.5*(hi-lo)*w



    elif floatmode == "mpfloat":

        efloat_mode="mpfloat"
        efloat = ef.EFloat(efloat_mode)

        T = efloat.fromArray(numpy.zeros([n,n]))

        for i in range(1,n):
            beta = .5 / efloat.sqrt(1.-(2*(i))**(-2.))   # 3-term recurrence coeffs

            # upper diagonal
            T[i,i-1] = beta

            # lower diagonal
            T[i-1,i] = beta

        if radau:
            # Eq. 2.4 of Gautschi, Gauss–Radau formulae for Jacobi and Laguerre weight functions, 2000
            if n == 1:
                T[-1,-1] = 1
            else:
                T[-1,-1] = 1 - 2*(n-1)**2 / (2*(n-1)*(2*(n-1)+1))

        D, V = efloat.eigh(T)                      # Eigenvalue decomposition

        N = len(D)

        i = numpy.argsort(D)         # Legendre points

        x = efloat.fromArray(numpy.zeros(N))
        w = efloat.fromArray(numpy.zeros(N))

        for i in range(len(D)):
            x[i] = D[i]

            w[i] = 2*V[0,i]**2

            x[i] = 0.5*((hi+lo) + (hi-lo)*x[i])
            w[i] = 0.5*(hi-lo)*w[i]

        x = numpy.array(x.tolist(), dtype=float)[:,0]
        w = numpy.array(w.tolist(), dtype=float)[:,0]

    else:
        raise Exception("Unsupported floatmode")

    return x, w


def lagrange(c):
    """Construct the Lagrange polynomials on the set of distinct abscissa
    c[].  The returned function can be evaluated at a vector of target
    locations, producing a Vandermonde-type matrix with columns of
    Lagrange polynomials and rows for each target point.
    """
    def fj(j, t):
        cc = numpy.concatenate((c[:j], c[j+1:]))
        return numpy.prod(t - cc) / numpy.prod(c[j] - cc)
    def f(targets):
        return numpy.array([ [fj(j,t) for j in range(len(c))] for t in targets])
    return f



def irk_collocation(c, golubwelsch=_golubwelsch):
    """Construct the Butcher table associated with the implicit Runge-Kutta
    collocation method on the points c.  We compute the necessary
    integrals using Gauss numerical quadrature.
    """
    s = len(c)
    A = numpy.empty((s,s), dtype=c.dtype)
    p = lagrange(c)
    x, w = golubwelsch(s, 0, 1)
    for i in range(s):
        A[i,:] = c[i]*w.dot(p(c[i]*x))
    b = w.dot(p(x))
    return A, b, c



def _tanhsinh(s, h):
    N = s

    k = numpy.array(range(-N//2, N//2+1))
    if floatmode == "float":

        x = numpy.tanh(0.5*numpy.pi*numpy.sinh(k*h))
        w = 0.5*h*numpy.pi*numpy.cosh(k*h)/numpy.power(numpy.cosh(0.5*numpy.pi*numpy.sinh(k*h)), 2.0)

        # To interval [0;1]
        x = (x+1.0)*0.5
        w *= 0.5

    elif floatmode == "mpfloat":

        efloat_mode="mpfloat"
        efloat = ef.EFloat(efloat_mode)
        import mpmath as mp

        x = efloat.fromArray(numpy.zeros([N]))
        w = efloat.fromArray(numpy.zeros([N]))

        for i in range(N):
            x[i,0] = mp.mp.tanh(0.5*mp.mp.pi*mp.mp.sinh(k[i]*h))
            w[i,0] = 0.5*h*mp.mp.pi*mp.mp.cosh(k[i]*h)/mp.mp.power(mp.mp.cosh(0.5*mp.mp.pi*mp.mp.sinh(k[i]*h)), 2.0)

        # To interval [0;1]
        x = (x+1.0)*0.5
        w *= 0.5

        x = numpy.array(x.tolist(), dtype=float)[:,0]
        w = numpy.array(w.tolist(), dtype=float)[:,0]

    return x, w



def tanhsinh(s, h=None):

    if h==None:
        # Larger h's lead to stronger clustering at the interval ends
        #h = numpy.pi*2.0*0.5/s

        # Tanh-Sinh High-Precision Quadrature, David H. Bailey, 19 Jan 2006
        # We use 2 instead of 12 since we don't need 1000 digits accuracy ;-)
        h = 1.0/(2.0**2)

    c, _ = _tanhsinh(s, h)

    return irk_collocation(c)



def gauss(s):
    c, _ = _golubwelsch(s, 0, 1)
    return irk_collocation(c)

def radau(s):
    c, _ = _golubwelsch(s, 0, 1, radau=True)
    return irk_collocation(c)

def chebyshev(s, left=False, right=False):
    l = numpy.arange(s)
    if not left:
        l += 1
    d = max(l)
    if not right:
        d += 1
    c = .5 - .5*numpy.cos(numpy.pi*l/d)
    return irk_collocation(c)


def butcher_diag(A, b):
    s = len(b)
    Lambda, V = numpy.linalg.eig(A)
    btilde = numpy.linalg.solve(V, numpy.ones(s))*(V.T.dot(b))
    return numpy.diag(Lambda), btilde, Lambda


def butcher2mso(A, b):
    M = A.copy()
    L = 0*M
    ell = numpy.linalg.solve(A, b)
    mu = 0*ell
    return M, mu, L, ell


def butcher2rexi(A, b, classic=False):
    A, b, c = butcher_diag(A, b)
    M, mu, L, ell = butcher2mso(A, b)
    if not numpy.isclose(ell.sum(), 1) and classic:
        raise RuntimeError('Cannot convert from MSO to classic REXI because sum(ell)=={}, but must be 1'.format(ell.sum()))
    alpha = -1/numpy.diag(M)
    beta = alpha * ell
    return alpha, beta
