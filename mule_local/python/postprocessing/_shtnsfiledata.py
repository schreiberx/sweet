#! /usr/bin/env python3
"""
DEPRECATED
DEPRECATED
DEPRECATED

Use SphereDataOperators

DEPRECATED
DEPRECATED
DEPRECATED
"""

import math, sys
import shtns
import numpy as np

class shtnsfiledata:
    """
    Adopted from https://bitbucket.org/nschaeff/shtns/src/master/examples/shallow_water.py
    """
    def __init__(
            self,
            rsphere = 1.0
    ):
        print("*"*80)
        print("DEPRECATED, don't use this, use SphereDataOperators instead!")
        print("*"*80)

        self.rsphere = rsphere


    def setup(self, file_info, anti_aliasing=False):
        import shtns
        import numpy as np

        if file_info['modes_m_max'] != file_info['modes_m_max']:
            raise Exception("Only num_lon == num_lat supported")

        ntrunc = file_info['modes_n_max']
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)

        nlons = (ntrunc + 1) * 2
        nlats = (ntrunc + 1)

        if anti_aliasing:
            if nlons & 1:
                raise Exception("Only even numbers of longitudinal coordinates allowed for anti-aliasing")
            if nlats & 1:
                raise Exception("Only even numbers of latitudinal coordinates allowed for anti-aliasing")

            print("Anti-aliasing:")
            print(" + old lon/lat: ", nlons, nlats)

            nlons += nlons//2
            nlats += nlats//2

            print(" + new lon/lat: ", nlons, nlats)

        if file_info['grid_type'] == 'GAUSSIAN':
            #self._shtns.set_grid(nlats,nlons,shtns.sht_gauss_fly|shtns.SHT_PHI_CONTIGUOUS, 1.e-10)
            self._shtns.set_grid(nlats, nlons, shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS, 0)

        elif file_info['grid_type'] == 'REGULAR':
            #self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS, 1.e-10)
            self._shtns.set_grid(nlats, nlons, shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS, 0)

        else:
            raise Exception("Grid type '"+file_info['grid_type']+"' not supported!")

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
        self.lap = self.lap/self.rsphere**2
        self.invlap = self.invlap*self.rsphere**2



    def phys2spec(self, data):
        return self._shtns.analys(data)

    def spec2phys(self, dataspec):
        return self._shtns.synth(dataspec)

    def vrtdiv2uv(self, vrtspec, divspec):
        return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)

    def uv2vrtdiv(self,u,v):
        vrtspec, divspec = self._shtns.analys(u, v)
        return self.lap*self.rsphere*vrtspec, self.lap*self.rsphere*divspec

    def getuv(self,divspec):
        vrtspec = np.zeros(divspec.shape, dtype=np.cdouble)
        u,v = self._shtns.synth(vrtspec,divspec)
        return u, v

    def rotateX90(self, i_field):
        return self._shtns.Xrotate90(i_field)

    def rotateY90(self, i_field):
        return self._shtns.Yrotate90(i_field)

    def rotateZ90(self, i_field, angle):
        return self._shtns.Zrotate(i_field, angle)

