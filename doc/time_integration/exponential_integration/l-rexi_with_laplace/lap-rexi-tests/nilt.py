import sys
import os
import warnings
import numpy as np 
import math

import sympy as sy
#from sympy import series, Symbol
#from sympy.functions import sin, cos, exp

class nilt:
    def __init__(self, contour="ellipse"):
        print("Contour: "+contour)
        self.contour=contour #circle or ellipse
        if contour=="ellipse":
            self.ellipse=True
            self.circle=False
        else:
            self.ellipse=False
            self.circle=True

    def set_parameters(self, gamma=2.0, delta=0.5):
        #Set parameters
        self.gamma = gamma
        if self.ellipse:
            self.delta = delta
            self.r = math.sqrt((1-self.delta/self.gamma)/(1+self.delta/self.gamma))
            self.d = self.r+1/self.r
        else:
            self.r = self.gamma
        return
    
    def set_quadrature(self, N, ntrunc = 0 ):
        self.N = N
        self.twopi = 2.0*math.pi
        self.theta = np.linspace(0,self.twopi, N, endpoint = False)
        
        self.expplustheta = np.exp(1j*self.theta)
        if self.ellipse:
            self.expminustheta = np.exp(-1j*self.theta)
            self.sn = (self.gamma * 1j / self.d) * (self.r*self.expplustheta + (1/self.r)*self.expminustheta )
            self.spn = (-self.gamma / self.d) * (self.r*self.expplustheta - (1/self.r)*self.expminustheta)
        else:
            self.sn = self.r*self.expplustheta
            self.spn = 1j*self.sn

        self.use_exp_trunc = False
        if ntrunc > 0:
            x = sy.Symbol('x')
            self.x=x
            exp_taylor =sy.functions.exp(x).series(x, 0, ntrunc).removeO()
            self.exp_trunc = sy.utilities.lambdify(x, exp_taylor,'numpy')
            self.use_exp_trunc = True

        return

    def apply_nilt(self, fhat, t):
        if self.use_exp_trunc:
            exp_tsn=self.exp_trunc(t*self.sn)
            #print(np.abs(exp_tsn-np.exp(t*self.sn)).sum())
        else:
            exp_tsn=np.exp(t*self.sn)
        
        self.nilt=exp_tsn*fhat(self.sn)*self.spn

        self.niltsum=(1/(self.N*1j))*self.nilt.sum()
        if self.ellipse:
            self.niltsum=-self.niltsum
        return self.niltsum


