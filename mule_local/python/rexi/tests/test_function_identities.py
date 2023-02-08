#! /usr/bin/env python3

import sys
import os
import numpy as np

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d+"/..")

from mule.rexi.Functions import *

sys.path.pop()

funs = Functions()
funs_mp = Functions(efloat_mode="mpfloat")

eps = 1e-10
eps_series = 1e-6
eps_direct = 1e-1

minmax = 2.0
numsamples=123*2


#
# Test A)
#
# phiNDirect == '''method'''
#

for method in ['phiNDirectFormula', 'phiNRec', 'phiNSeries', 'phiN']:
#for method in []:
    print("*"*80)
    print("* Testing "+method)
    print("*"*80)
    for n in range(0, 6):
        print("* Phi "+str(n))

        maxerr = 0
        for xr in np.linspace(-minmax, minmax, numsamples):
        #if True:
            xr = 0

            for xi in np.linspace(-minmax, minmax, numsamples):
                x = xr + xi*1j

                if method=="phiNSeries":
                    if np.abs(x) < eps_series:
                        continue

                    # Compuare Series formulation with multi-precision implementation
                    val = funs_mp.phiNDirect(n, x)
                    val = complex(val)
                    valf = getattr(funs, method)(n, x)
                else:
                    if np.abs(x) < eps_direct and method in ['phiNDirectFormula', 'phiNRec']:
                        continue

                    val = funs.phiNSeries(n, x)
                    val = complex(val)
                    valf = getattr(funs, method)(n, x)

                err = funs.efloat.abs(val-valf)

                maxerr = max(err, maxerr)
                #eps_threshold = eps*funs.efloat.pow(10, n*2)
                eps_threshold = eps #*funs.efloat.pow(10, n*2)

                if err > eps_threshold:
                    print(" + x: "+str(x))
                    print(" + val: "+str(val))
                    print(" + valf: "+str(valf))
                    print(" + err: "+str(err))
                    print(" + err_threshold: "+str(eps_threshold))
                    raise Exception("Error too high!")

        print(" + maxerr: "+str(maxerr))


#
# Test B)
#
# upsNDirect == '''method'''
#

for method in ['upsNSeries', 'upsN']:
    print("*"*80)
    print("* Testing "+method)
    print("*"*80)
    for n in range(1, 4):
        print("* Ups "+str(n))

        maxerr = 0
        for xr in np.linspace(-minmax, minmax, numsamples):
            xr = 0

            for xi in np.linspace(-minmax, minmax, numsamples):
                x = xr + xi*1j


                if method=="upsNSeries":
                    if np.abs(x) < 1e-2:
                        continue

                    val = funs_mp.upsNDirect(n, x)
                    val = complex(val)
                    valf = getattr(funs, method)(n, x)

                else:
                    if np.abs(x) < eps_direct:
                        continue

                    val = funs.upsNSeries(n, x)
                    val = complex(val)
                    valf = getattr(funs, method)(n, x)

                err = funs.efloat.abs(val-valf)

                maxerr = max(err, maxerr)

                if err > eps:
                    print(" + x: "+str(x))
                    print(" + val: "+str(val))
                    print(" + valf: "+str(valf))
                    print(" + err: "+str(err))
                    raise Exception("Error too high!")

        print(" + maxerr: "+str(maxerr))


