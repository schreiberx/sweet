#! /usr/bin/env python3

import shtns


class shtnsfiledata:
    #
    # Adopted from https://bitbucket.org/nschaeff/shtns/src/master/examples/shallow_water.py
    #
    def __init__(
            self,
            rsphere = 1.0
    ):
        self.rsphere = rsphere


    def setup(self, file_info):
        import shtns
        import numpy as np

        if file_info['modes_m_max'] != file_info['modes_m_max']:
            raise Exception("Only num_lon == num_lat supported")

        ntrunc = file_info['modes_n_max']
        self._shtns = shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)

        nlons = (ntrunc+1)*2
        nlats = (ntrunc+1)

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
        self.lap = -self.degree*(self.degree+1.0).astype(np.complex)
        self.invlap = np.zeros(self.lap.shape, self.lap.dtype)
        self.invlap[1:] = 1./self.lap[1:]
        self.lap = self.lap/self.rsphere**2
        self.invlap = self.invlap*self.rsphere**2


    def phys2spec(self, data):
        return self._shtns.analys(data)

    def spec2phys(self, dataspec):
        return self._shtns.synth(dataspec)

    def vortdiv2uv(self, vrtspec, divspec):
        return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)

    def uv2vortdiv(self,u,v):
        vrtspec, divspec = self._shtns.analys(u, v)
        return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec

    def getuv(self,divspec):
        vrtspec = np.zeros(divspec.shape, dtype=np.complex)
        u,v = self._shtns.synth(vrtspec,divspec)
        return u, v

