#! /usr/bin/env python3

#
# Author: Martin Schreiber
# Email: M.Schreiber@exeter.ac.uk
# Date: 2017-06-10
#

import sys


default_floatmode = 'float'
default_mpfloat_accuracy = 60


class EFloat:
    __shared_state = {}

    def __init__(self, floatmode = None, accuracy = None):
        global default_floatmode, default_mpfloat_accuracy

        if floatmode == None:
            self.floatmode = default_floatmode
        else:
            self.floatmode = floatmode
            default_floatmode = floatmode

        if accuracy == None:
            self.accuracy = default_mpfloat_accuracy
        else:
            self.accuracy = accuracy
            default_mpfloat_accuracy = accuracy

        if self.floatmode in ['float']:
            import numpy as np
            import scipy as sp

            def to(value):
                if isinstance(value, list):
                    return [float(x) for x in value]
                return float(value)
            self.to = to

            def toCplx(value):
                if isinstance(value, list):
                    return [complex(x) for x in value]
                return complex(value)
            self.toCplx = toCplx

            self.linspace = lambda xmin, xmax, M: np.linspace(xmin, xmax, M).tolist()
            self.arange = lambda xmin, xmax, M: np.arange(xmin, xmax, M).tolist()

            self.i = 1j

            self.exp = lambda x : np.exp(x)
            self.pow = lambda x, y : np.power(x, y)
            self.sqrt = lambda x : np.sqrt(x)

            self.pi = np.pi
            self.pi2 = 2.0*np.pi

            def max_(a, b):
                return max(a, b)
            self.max = max_

            def abs_(a):
                return abs(a)
            self.abs = abs_

            def re(value):
                if isinstance(value, list):
                    return [x.real for x in value]
                return value.real
            self.re = re

            def im(value):
                if isinstance(value, list):
                    return [x.imag for x in value]
                return value.imag
            self.im = im

            self.conj = lambda value: value.real - 1j*value.imag

            def strx(value):
                return str(value)
            self.str = strx

            def floatToStr(value):
                return f'{value:.20f}'
            self.floatToStr = floatToStr

            self.zeroMatrixReal = lambda sx, sy: np.zeros((sx, sy), dtype='float64')

            self.zeroMatrixReal = lambda sx, sy: np.zeros((sx, sy), dtype='float64')
            self.zeroMatrixComplex = lambda sx, sy: np.zeros((sx, sy), dtype='complex128')

            def matvecmul(M, b):
                return M.dot(b)
            self.matvecmul = matvecmul

            def eigh(M):
                return np.eigh(M)
            self.eigh = eigh


            def transpose(M):
                return M.transpose()
            self.transpose = transpose

            def transpose_conj(M):
                return M.transpose_conj()
            self.transpose_conj = transpose_conj

            def matmatmul(Ma, Mb):
                return np.matmul(Ma, Mb)
            self.matmatmul = matmatmul

            def testHermitian(M):
                print(M.shape)
                T = M * np.conj(M.transpose()) - np.diag([1.0 for i in range(max(M.shape[0], M.shape[1]))])
                print(T)

                maxvalre = 0
                maxvalim = 0
                for j in range(T.shape[0]):
                    for i in range(T.shape[1]):
                        x = T[j,i]
                        maxvalre = max(abs(x.real), maxvalre)
                        maxvalim = max(abs(x.imag), maxvalim)

                print("maxvalre: "+str(float(maxvalre)))
                print("maxvalim: "+str(float(maxvalim)))

            self.testHermitian = testHermitian


            def psolve(M, b):
                if False:
                #if True:
                    # Only to test SVD
                    U, S, V = np.linalg.svd(M, compute_uv=True, full_matrices=True)

                    print("V")
                    testHermitian(V)

                    print("U")
                    testHermitian(U)

                    sys.exit(1)

                else:
                    Minv = np.linalg.pinv(M)
                ws = Minv.dot(b)
                return ws
            self.psolve = psolve 

            def solve(M, b):
                ws = np.linalg.solve(M, b)
                return ws
            self.solve = solve 

            def lu_solve_lsq(M, b):
                print("LU decomposition for least squared problem not available with mode 'float'")
                return

            self.lu_solve_lsq = lu_solve_lsq

            def qr_solve_lsq(M, b):
                Q, R = np.linalg.qr(M)
                Vt = np.linalg.inv(R)*(np.conj(Q.transpose())*np.matrix(b).transpose())

                Vt = Vt.tolist()
                V = np.zeros(M.shape[1], dtype=M.dtype)
                for i in range(M.shape[1]):
                    V[i] = Vt[i][0]

                return V

            self.qr_solve_lsq = qr_solve_lsq

            def tolist(x):
                return [i for i in x]
            self.tolist = tolist

            def norm2(v):
                return np.linalg.norm(v, 2)
            self.norm2 = norm2


            def tovector(x):
                return np.matrix(x).transpose()
            self.tovector = tovector

            self.ceil = lambda value: np.ceil(value)
            self.floor = lambda value: np.floor(value)

            self.formatFloat = lambda value, digits: ("%0."+str(digits)+"e") % value


        elif self.floatmode == 'mpfloat':
            import mpmath as mp

            mp.mp.dps = self.accuracy

            def fromArray(M):
                return mp.mp.matrix(M)
            self.fromArray = fromArray


            def to(value):
                if isinstance(value, list):
                    return [mp.mp.mpf(x) for x in value]
                return mp.mp.mpf(value)
            self.to = to

            def toCplx(value):
                if isinstance(value, list):
                    return [mp.mp.mpc(x) for x in value]
                return mp.mp.mpc(value)
            self.toCplx = toCplx

            self.linspace = lambda xmin, xmax, M: mp.mp.linspace(xmin, xmax, M)
            self.arange = lambda xmin, xmax, M: mp.mp.arange(xmin, xmax, M)

            self.i = mp.mp.j

            self.exp = lambda x : mp.mp.exp(x)
            self.pow = lambda x, y : mp.mp.power(x, y)
            self.sqrt = lambda x : mp.mp.sqrt(x)

            self.pi = mp.mp.pi
            self.pi2 = mp.mp.pi*self.to(2.0)

            def max_(a, b):
                return mp.mp.max(a, b)
            self.max = max_

            def abs_(a):
                return abs(a)
            self.abs = abs_

            def re(value):
                if isinstance(value, list):
                    return [mp.mp.re(x) for x in value]
                return mp.mp.re(value)
            self.re = re

            def im(value):
                if isinstance(value, list):
                    return [mp.mp.im(x) for x in value]
                return mp.mp.im(value)
            self.im = im

            self.conj = lambda value: mp.mp.re(value) - mp.mp.j*mp.mp.im(value)

            self.str = lambda var: str(var)

            def floatToStr(value):
                return str(value)
            self.floatToStr = floatToStr

            self.zeroMatrixReal = lambda sx, sy: mp.mp.zeros(sx, sy)
            self.zeroMatrixComplex = lambda sx, sy: mp.mp.zeros(sx, sy)


            def printM(M):
                for j in range(M.rows):
                    for i in range(M.cols):
                        v = M[j,i]
                        if abs(v) < 1e-10:
                            v = 0.0
                        print(str(v)+" ", end="")
                    print("")

            def testHermitian(M):
                T = M * (M.transpose_conj()) - mp.diag([1.0 for i in range(max(M.rows, M.cols))])
                #printM(T)

                maxvalre = 0
                maxvalim = 0
                for j in range(T.rows):
                    for i in range(T.cols):
                        x = T[j,i]
                        maxvalre = max(abs(x.real), maxvalre)
                        maxvalim = max(abs(x.imag), maxvalim)

                print("maxvalre: "+str(float(maxvalre)))
                print("maxvalim: "+str(float(maxvalim)))
            self.testHermitian = testHermitian


            def psolve(M, b):
                U, S, V = mp.svd(M)

                if False:
                    print("V")
                    testHermitian(V)
                    print("U")
                    testHermitian(U)
                    sys.exit(1)

                Mpinv = V.transpose_conj() * (mp.diag(S)**(-1)) * U.transpose_conj()
                ws = Mpinv*mp.matrix(b)
                return [ws[i] for i in range(len(ws))]

            self.psolve = psolve

            def transpose(M):
                return M.transpose()
            self.transpose = transpose

            def transpose_conj(M):
                return M.transpose_conj()
            self.transpose_conj = transpose_conj

            def matmatmul(Ma, Mb):
                return Ma*Mb
            self.matmatmul = matmatmul

            def matvecmul(M, b):
                return M*mp.mp.matrix(b)
            self.matvecmul = matvecmul

            def eigh(M):
                return mp.mp.eigh(M)
            self.eigh = eigh

            def lu_solve_lsq(M, b):
                x = mp.lu_solve(M, b)
                return x

            self.lu_solve_lsq = lu_solve_lsq


            def solve(M, b):
                return M**(-1)*b
            self.solve = solve

            def qr_solve_lsq(M, b):
                Vt = mp.qr_solve(M, b)

                V = mp.zeros(M.cols,1)
                for i in range(M.cols):
                    V[i] = Vt[0][i]

                return V


            self.qr_solve_lsq = qr_solve_lsq

            def tolist(x):
                return [i for i in x]
            self.tolist = tolist

            self.ceil = lambda value: mp.ceil(value)
            self.floor = lambda value: mp.floor(value)

            self.formatFloat = lambda value, digits: ("%0."+str(digits)+"e") % float(value)

        else:
            raise Exception("Unknown floatmode "+str(self.floatmode))

        self.cplx = lambda re, im: self.to(re) + self.i*self.to(im)

