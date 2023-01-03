#! /usr/bin/env python3

import numpy as np
import shtns
import sys


class Spharmt(object):
    def __init__(self,nlons,nlats,ntrunc,rsphere,gridtype='gaussian'):
        """initialize
        nlons:  number of longitudes
        nlats:  number of latitudes
        """
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
        if gridtype == 'gaussian':
            self._shtns.set_grid(nlats,nlons,shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,0)
        elif gridtype == 'regular':
            self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,0)
        self.lats = np.arcsin(self._shtns.cos_theta)
        self.lons = (2.*np.pi/nlons)*np.arange(nlons)
        
    def getvrtdivspec(self,u,v):
        vrtspec, divspec = self._shtns.analys(u, v)
        return None, None


print("This is a segfault test!")

#
# Comment this out to avoid segfault!
#
print("Setting up non-aliased formulation")
y = Spharmt(2048,1024,1023,rsphere=1,gridtype='gaussian')

# Setup aliased version
print("Setting up aliased formulation")
ntrunc = 1023
nlons = 3072
nlats = 1536
x = Spharmt(nlons,nlats,ntrunc,rsphere=1,gridtype='gaussian')
lons,lats = np.meshgrid(x.lons, x.lats)

ug = np.zeros((nlats,nlons),np.double)
vg = np.zeros((nlats,nlons),np.double)

print("Running transformation")
# Do transformation (should segfault)
vrtspec, divspec =  x.getvrtdivspec(ug,vg)

print("Everything OK")
