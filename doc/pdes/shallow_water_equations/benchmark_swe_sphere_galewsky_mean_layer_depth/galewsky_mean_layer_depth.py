#! /usr/bin/env python3

#
# Compute mean layer depth for Galewsky benchmark
#
# See Eq. (3) and comment on computing h0 on p. 431, left column top
#

import numpy as np
import shtns
import sys


class Spharmt(object):
    """
    wrapper class for commonly used spectral transform operations in
    atmospheric models.  Provides an interface to shtns compatible
    with pyspharm (pyspharm.googlecode.com).
    """
    def __init__(self,nlons,nlats,ntrunc,rsphere,gridtype='gaussian',verbose=False):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes"""
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
        if gridtype == 'gaussian':
            self._shtns.set_grid(nlats,nlons,shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,0)
        elif gridtype == 'regular':
            self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,0)
        self.lats = np.arcsin(self._shtns.cos_theta)
        self.lons = (2.*np.pi/nlons)*np.arange(nlons)
        self.nlons = nlons
        self.nlats = nlats
        self.ntrunc = ntrunc
        self.nlm = self._shtns.nlm
        self.degree = self._shtns.l
        self.lap = -self.degree*(self.degree+1.0).astype(np.cdouble)
        self.invlap = np.zeros(self.lap.shape, self.lap.dtype)
        self.invlap[1:] = 1./self.lap[1:]
        self.rsphere = rsphere
        self.lap = self.lap/rsphere**2
        self.invlap = self.invlap*rsphere**2
        
        if verbose:
            print("N: "+str(self.nlons)+", "+str(self.nlats))
            print("Mtrunc: "+str(self.ntrunc))
            print("Nlm: "+str(self.nlm))

    def grdtospec(self,data):
        """compute spectral coefficients from gridded data"""
        return self._shtns.analys(data)
    def spectogrd(self,dataspec):
        """compute gridded data from spectral coefficients"""
        return self._shtns.synth(dataspec)
    def getuv(self,vrtspec,divspec):
        """compute wind vector from spectral coeffs of vorticity and divergence"""
        return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)
    def getvrtdivspec(self,u,v):
        """compute spectral coeffs of vorticity and divergence from wind vector"""
        vrtspec, divspec = self._shtns.analys(u, v)
        return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec






#
# Constants from the paper
#

a = 6.37122e6   # earth radius
omega = 7.292e-5    # rotation rate
g = 9.80616  # gravity
hbar = 10.e3    # resting depth
umax = 80.  # jet speed
phi0 = np.pi/7.
phi1 = 0.5*np.pi - phi0
phi2 = 0.25*np.pi
en = np.exp(-4.0/(phi1-phi0)**2)
alpha = 1./3.; beta = 1./15.
hamp = 120.



#
# Velocity profile computed at phi
#
# See Eq. (2)
#
def u_init(phi):
    if 0:
        # Numerically not that stable method
        u_ = (umax/en)*np.exp(1./((phi-phi0)*(phi-phi1)))
        return np.where(np.logical_and(phi < phi1, phi > phi0), u_, np.zeros_like(phi))
    else:
        # Much better, don't get close to 1/0
        e = 1e-8
        cond = np.logical_and(phi < phi1-e, phi > phi0+e)
        denom = np.where(cond, (phi-phi0)*(phi-phi1), np.ones_like(phi))
        u_ = (umax/en)*np.exp(1./denom)
        return np.where(cond, u_, np.zeros_like(phi))



eps = 1e-16


# Return circumfence of circle at latituce
def metric_term(lat):
    return np.cos(lat)*2.0*np.pi


#
# Integral in Eq. 3
#
def int_f(lat):
    f = 2.*omega*np.sin(lat)
    return a*u_init(lat)*(f+np.tan(lat)/a*u_init(lat))


import scipy.integrate as integrate

#
# Quadrature of integral int_f in Eq. 3
#
def quad_f(lat):

    if not isinstance(lat, np.ndarray):
        lat = np.array([lat])

    ret = np.zeros_like(lat)
    for i in range(len(lat)):
        (quad, _) = integrate.quad(lambda x: int_f(x), 0, lat[i], epsrel=eps)
        ret[i] = quad

    return ret




print("*"*80)
print("Min/max of quadratured function / g")
print(" + quad_f/g min: "+str(quad_f(-np.pi*0.5)/g))
print(" + quad_f/g max: "+str(quad_f(np.pi*0.5)/g))


#
# Compute mean layer depth h0
# (p. 430, left upper column, Galewsky et al. paper)
#

# Compute over entire sphere's surface
(quad_total, _) = integrate.quad(lambda x: quad_f(x)*metric_term(x), -0.5*np.pi, 0.5*np.pi, epsrel=eps)
print("*"*80)
print("quad of quad_f \in [-pi/2;pi/2]: "+str(quad_total))

# Normalize across quadrature region
gh_avg = quad_total/(np.pi*4)

h_avg = gh_avg/g
print(" + h_avg: "+str(h_avg))


# We use a + sign here, since the integral is subtracted
h0 = 10000 + h_avg

print("*"*80)
print("*"*80)
print("*"*80)
print("*"*80)
print("* THIS IS THE GALEWSKY BENCHMARK H0 CONSTANT!!!!!!!!!!!")
print("* h0: "+str(h0))
print("*"*80)
print("*"*80)
print("*"*80)
print("*"*80)
print("*"*80)



#
# Test routine
#

for n in [2**i for i in range(5, 12)]:

    print("*"*80)
    print("* Testing with resolution "+str(n))

    nlons = n
    nlats = n//3*2
    ntrunc = int(nlons/3)  # spectral truncation (for alias-free computations)

    x = Spharmt(nlons,nlats,ntrunc,a,gridtype='gaussian',verbose=False)

    data = h0*g - quad_f(x.lats)

    print("*"*80)
    print(" + min(data): "+str(np.min(data)))
    print(" + max(data): "+str(np.max(data)))


    # Use lats from negative to positive numbers
    lats = np.flip(x.lats)
    areas = np.zeros(lats.shape[0])

    dlat = (lats[1]+lats[0])*0.5 + np.pi*0.5
    assert dlat > 0 and dlat < 0.5
    areas[0] = dlat*metric_term(lats[0])

    for i in range(1, len(areas)-1):
        dlat = (lats[i+1]-lats[i-1])*0.5
        assert dlat > 0 and dlat < 0.5

        areas[i] = dlat*metric_term(lats[i])

    dlat = np.pi*0.5 - (lats[-2]+lats[-1])*0.5
    assert dlat > 0 and dlat < 0.5

    areas[-1] = dlat*metric_term(lats[-1])

    #areas /= 2.0

    # Normalize along latitude
    areas /= 4.0*np.pi

    avg_gh = np.sum(areas*data)
    avg_h = avg_gh/g
    print(" + numerical avg h: "+str(avg_h))
    print("")
