#! /usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable


input_file = sys.argv[1]

f = open(input_file, 'rb')
content = f.read()
f.close()

numbers = re.findall("\d+", input_file)
time = int(numbers[-2])/24

file_info = {}

acc = ""
line_nr = 0
fin_detected = False

i = 0
for i in range(0, len(content)):
    c = content[i]

    if c == ord('\n'):
        if line_nr == 0:
            if acc == "SWEET":
                print("SWEET header detected")
            else:
                raise Exception("Header not detected, stopping")

        else:
            s = "DATA_TYPE"
            if acc[0:len(s)] == s:
                file_info['data_type'] = acc[len(s)+1:]

            s = "MODES_N_MAX"
            if acc[0:len(s)] == s:
                file_info['modes_n_max'] = int(acc[len(s)+1:])

            s = "MODES_M_MAX"
            if acc[0:len(s)] == s:
                file_info['modes_m_max'] = int(acc[len(s)+1:])

            s = "NUM_ELEMENTS"
            if acc[0:len(s)] == s:
                file_info['num_elements'] = int(acc[len(s)+1:])

            s = "GRID_TYPE"
            if acc[0:len(s)] == s:
                file_info['grid_type'] = acc[len(s)+1:]

            s = "FIN"
            if acc[0:len(s)] == s:
                print("FIN detected")
                fin_detected = True
                i += 1
                break

        acc = ""
        line_nr += 1

    else:
        acc += chr(c)

if not fin_detected:
    raise Exception("FIN not detected in file")

if file_info['data_type'] != "SH_DATA":
    raise Exception("Wrong data type "+file_info['data_type']+" in binary file")


print("*"*80)
for key, value in file_info.items():
    print(str(key)+": "+str(value))
print("*"*80)


data = content[i:]
print("BINARY DATA LEN: "+str(len(data)))


#
# SHTNS STUFF HERE
#

class SHTNS_data:
    #
    # Adopted from https://bitbucket.org/nschaeff/shtns/src/master/examples/shallow_water.py
    #
    def __init__(
            self,
            rsphere = 6.37122e6
    ):
        self.rsphere = rsphere

    def setup(self, file_info, data):
        import shtns
        import numpy as np

        if file_info['modes_m_max'] != file_info['modes_n_max']:
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

        sh_data = np.frombuffer(data, dtype=np.complex128)

        return sh_data


    def phys2spec(self, data):
        return self._shtns.analys(data)

    def spec2phys(self, dataspec):
        return self._shtns.synth(dataspec)

    def vortdiv2uv(self, vrtspec, divspec):
        return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)

    def uv2vortdiv(self,u,v):
        vrtspec, divspec = self._shtns.analys(u, v)
        return self.lap*self.rsphere*vrtspec, self.lap*self.rsphere*divspec

    def getuv(self,divspec):
        vrtspec = np.zeros(divspec.shape, dtype=np.complex)
        u,v = self._shtns.synth(vrtspec,divspec)
        return u, v


s = SHTNS_data()
data_spec = s.setup(file_info, data)

data_phys = s.spec2phys(data_spec)

fig = plt.figure(figsize=(8, 4))
ax = plt.gca()
# dimensionless PV
lons1d = (180./np.pi)*s.lons-180.
lats1d = (180./np.pi)*s.lats
lmax = np.max(np.max(data_phys))
lmin = np.min(np.min(data_phys))
ld = (lmax-lmin)/50
levs = np.arange(lmin, lmax, ld)
print(levs)

#cs = plt.contourf(lons1d, lats1d, data_phys, levs, extend="both", aspect='auto', cmap="nipy_spectral")
if 'vrt' in input_file:
    levs = np.arange(-7e-6, 7e-6, 1e-6)
    cs = plt.contourf(lons1d, lats1d, data_phys, levs, extend="both", cmap="bwr") #, linewidth=0.01)
else:
    levs = np.arange(lmin, lmax, ld)
    cs = plt.contourf(lons1d, lats1d, data_phys, levs, extend="both", cmap="nipy_spectral") #, linewidth=0.01)

for c in cs.collections:
    c.set_edgecolor("face")

plt.grid(c='gray', ls='--', alpha=0.3)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.xticks(np.arange(-180, 180.1, 60))
plt.yticks(np.arange(-90, 90.1, 30))
plt.title(r'$T = %3.1f \, \mathrm{days}, \, \alpha=20.0$'%(time))


divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)
cb = plt.colorbar(cs, orientation="vertical", cax=cax) #, aspect=20, shrink=0.4, pad=0.2, cax=cax)
cb.set_label("Vorticity")

#plt.axis("equal")
#plt.axis("tight")


plt.tight_layout()

outputfile = input_file.replace('.sweet', '.png')

if input_file == outputfile:
    raise Exception("Input file didn't end with .sweet")

print("Writing to "+str(outputfile))
plt.savefig(outputfile)

