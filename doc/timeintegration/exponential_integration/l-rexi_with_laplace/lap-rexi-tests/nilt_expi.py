
import sys
import os
import warnings
import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

from nilt import nilt

#function reconstrution
def nilt_recon(alpha, ntrunc=0):

    #Test function
    def fhat(s):
        return 1/(s-1j*alpha)

    def f(t):
        return np.exp(1j*t*alpha)

    nend = 512
    nini = 32
    n_points=np.arange(nini, nend, 4)
    print(n_points)

    ndt = 40
    dtini=0.0
    dtend=20.0
    dtstep = (dtend-dtini)/ndt
    dt_points=np.arange(dtini, dtend, dtstep)
    print(dt_points)

    error_circ=np.zeros((len(dt_points),len(n_points)))
    error_ellip=np.zeros((len(dt_points),len(n_points)))

    print("dt, N Circ_Error Ellip_Error")
    for i, N in enumerate(n_points):
        if ntrunc < 0:
            M=N
        else:
            M=ntrunc
        circ.set_quadrature(N, ntrunc = M)
        ellip.set_quadrature(N, ntrunc = M)

        for j, dt in enumerate(dt_points):
            nilt_circ=circ.apply_nilt(fhat, dt)
            nilt_ellip=ellip.apply_nilt(fhat,dt)
            error_c=np.abs(np.abs(nilt_circ-f(dt)))
            error_e=np.abs(np.abs(nilt_ellip-f(dt)))    
            if error_e > 100:
                error_e = np.nan
            if error_c > 100:
                error_c = np.nan
            if error_c < 1e-15:
                error_c = 1e-16
            if error_e < 1e-15:
                error_e = 1e-16
            error_circ[j, i]=error_c
            error_ellip[j, i]=error_e
            print(dt, N, error_circ[j, i], error_ellip[j, i])
    print()
    return n_points, dt_points, error_circ, error_ellip

ellip = nilt("ellipse")
ellip.set_parameters(10.0, 2.0)

circ =  nilt("circle")
circ.set_parameters(10.0)

circ.set_quadrature(32)
ellip.set_quadrature(32)

fig1, ax = plt.subplots()
ax.scatter(ellip.sn.real,ellip.sn.imag, label="Ellipse", marker=".")
ax.scatter(circ.sn.real,circ.sn.imag, label="Circle", marker=".")
ax.legend(loc="best")
plt.savefig("ellip_circ.png")

fig3, axes = plt.subplots(4,2, figsize=(10,15))
plt.suptitle(r'$e^{i\alpha}$'+" reconstruction "+r'$(\gamma=$'+str(ellip.gamma)+", "+r'$\delta=$'+str(ellip.delta)+")")
levels = np.geomspace(1e-16, 1, num=17)
cmap=cm.viridis

alphalist=[0.0, 9.0, 9.9, 9.99]
ntrunc=0
for i, alpha in enumerate(alphalist):
    n_points, dt_points, error_circ, error_ellip = nilt_recon(alpha, ntrunc=ntrunc)
    X, Y = np.meshgrid(n_points, dt_points)
    cs = axes[i,0].contourf(X, Y, error_circ, levels, locator=ticker.LogLocator(), cmap=cmap)
    axes[i,0].set(xlabel="Number of quadrature points", 
        ylabel="dt", title="Circle "+r'$\alpha=$'+str(alpha))    
    cs = axes[i,1].contourf(X, Y, error_ellip, levels, locator=ticker.LogLocator(), cmap=cmap)
    axes[i,1].set(xlabel="Number of quadrature points", 
        ylabel="dt", title="Ellipse "+r'$\alpha=$'+str(alpha))

fig3.tight_layout(pad=4.0)
fig3.subplots_adjust(right=0.8)
fig3.subplots_adjust(top=0.9)
cbar_ax = fig3.add_axes([0.85, 0.15, 0.05, 0.7])
cbar=fig3.colorbar(cs, cax=cbar_ax)
cbar.set_label('Error in NILT')

if ntrunc > 0:
    plt.savefig("expi_recon_trunc"+str(ntrunc)+".png")
elif ntrunc < 0:
    plt.savefig("expi_recon_truncN.png")
else:
    plt.savefig("expi_recon.png")

plt.show()
    


