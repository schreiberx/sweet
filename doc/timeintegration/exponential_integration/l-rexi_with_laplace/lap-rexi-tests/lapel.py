
import sys
import os
import warnings
import numpy as np 
import math
import matplotlib.pyplot as plt

from nilt import *
from nilt_constant import *


#Constant function reconstrution
def constant_recon(dt, trunc=0):

    #Test function
    def fhat(s):
        return 1/s

    nend = 256
    nini = 32
    nstore =  int((nend-nini)/4)
    error_circ=np.zeros(nstore)
    error_ellip=np.zeros(nstore)
    quad_points=np.zeros(nstore)
    print("dt, N Circ_Error Ellip_Error")
    for i, N in enumerate(range(nini, nend, 4)):
        quad_points[i]=N
        circ.set_quadrature(N)
        ellip.set_quadrature(N)
        nilt_circ=circ.apply_nilt(fhat, dt, trunc)
        nilt_ellip=ellip.apply_nilt(fhat,dt, trunc)
        error_circ[i]=np.abs(np.abs(nilt_circ)-1.0)
        error_ellip[i]=np.abs(np.abs(nilt_ellip)-1.0)
        print(dt, N,error_circ[i], error_ellip[i])
    print()
    return quad_points, error_circ, error_ellip

ellip = nilt("ellipse")
ellip.set_parameters(10.0, 0.5)

circ =  nilt("circle")
circ.set_parameters(10.0)

circ.set_quadrature(32)
ellip.set_quadrature(32)

fig1, ax = plt.subplots()
ax.scatter(ellip.sn.real,ellip.sn.imag, label="Ellipse", marker=".")
ax.scatter(circ.sn.real,circ.sn.imag, label="Circle", marker=".")
ax.legend(loc="best")
plt.savefig("ellip_circ.png")

fig2, axes = plt.subplots(3,3, constrained_layout=True, figsize=(10,15))
plt.suptitle("Constant reconstruction")

dtlist = [0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0]

for  i, ax in enumerate(axes.reshape(-1)):
    dt=dtlist[i]
    quad_points, error_circ, error_ellip = constant_recon(dt, trunc=1)
    ax.plot(quad_points,error_circ, label="Circle")
    ax.plot(quad_points,error_ellip, label="Ellipse")
    ax.set(xlabel="Number of quadrature points", 
        ylabel="Error", yscale='log', title="dt"+str(dt))
    ax.legend(loc="best")

plt.savefig("const_recon_trunc.png")


plt.show()
    


