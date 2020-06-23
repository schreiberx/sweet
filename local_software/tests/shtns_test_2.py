#! /usr/bin/env python3

import sys
import numpy as np
import shtns


class shtnsfiledata:

    def setup(self, file_info, anti_aliasing):
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
                raise Exception("Only even numbers of longitudes coordinates allowed for anti-aliasing")
            if nlats & 1:
                raise Exception("Only even numbers of latitudinal coordinates allowed for anti-aliasing")

            print("Anti-aliasing:")
            print(" + old lon/lat: ", nlons, nlats)

            nlons += nlons//2
            nlats += nlats//2

            print(" + new lon/lat: ", nlons, nlats)

        self._shtns.set_grid(nlats, nlons, shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS, 0)

        self.lats = np.arcsin(self._shtns.cos_theta)
        self.lons = (2.*np.pi/nlons)*np.arange(nlons)



    def phys2spec(self, data):
        return self._shtns.analys(data)

    def spec2phys(self, dataspec):
        return self._shtns.synth(dataspec)



class postprocessing_swe:
    def run(self, M=128):
        file_info = {
                'modes_m_max': M-1,
                'modes_n_max': M-1,
                }

        #
        # Swap this between 0 and 1
        #
        if 1:
            #
            # Bug reproducer
            #

            # Setup transformations without anti-aliasing (lower physical resolution)
            self.sh = shtnsfiledata()
            self.sh.setup(file_info, anti_aliasing=False)


            # Setup transformations *with* anti-aliasing
            self.sh_aa = shtnsfiledata()
            self.sh_aa.setup(file_info, anti_aliasing=True)

        else:
            # Setup transformations *with* anti-aliasing
            self.sh_aa = shtnsfiledata()
            self.sh_aa.setup(file_info, anti_aliasing=True)

            # Setup transformations without anti-aliasing (lower physical resolution)
            self.sh = shtnsfiledata()
            self.sh.setup(file_info, anti_aliasing=False)


        # Load fields
        nlats = M
        nlons = M*2
        self.phi_pert_phys = np.zeros((nlats, nlons))

        phi0 = np.pi/7.
        phi1 = 0.5*np.pi - phi0
        phi2 = 0.25*np.pi
        alpha = 1./30.; beta = 1./150.
        hamp = 120.

        lons, lats = np.meshgrid(self.sh.lons, self.sh.lats)
        self.phi_pert_phys += hamp*np.cos(lats)*np.exp(-((lons-np.pi)/alpha)**2)*np.exp(-((phi2-lats)/beta)**2)

        self.phi_pert_spec = self.sh.phys2spec(self.phi_pert_phys)


        ke_phys_data = self.sh.spec2phys(self.phi_pert_spec)
        ke_spec_data = self.sh.phys2spec(ke_phys_data)
        err = np.max(np.abs(ke_spec_data - self.phi_pert_spec))
        print("Error", err)
        assert err < 1e-10

        ke_phys_data = self.sh_aa.spec2phys(self.phi_pert_spec)
        ke_spec_data = self.sh_aa.phys2spec(ke_phys_data)
        err = np.max(np.abs(ke_spec_data - self.phi_pert_spec))
        print("Error", err)
        assert err < 1e-10


print("Accuracy test")

for M in [64, 128, 256, 512, 1024, 2048]:
    s = postprocessing_swe()
    s.run()

print("Test successful")



