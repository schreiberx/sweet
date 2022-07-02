#! /usr/bin/env python3

import sys
import numpy as np

from mule.postprocessing.JobData import JobData
from mule_local.postprocessing.SphereDataSpectral import SphereDataSpectral
import mule_local.postprocessing.SphereDataOperators as SphereDataOperators
from postprocessing_spectrum_lib import postprocessing_spectrum_lib

debug_active = False

print(" ".join(sys.argv))

if len(sys.argv) < 5:
    print("Usage: "+sys.argv[0]+" 'jobdirectory' 'phidatafile.sweet' 'vrtdatafile.sweet' 'divdatafile.sweet'")
    raise Exception("Not enough arguments")

s = postprocessing_spectrum_lib()
s.setup(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

s.plot_kinetic_energy_spectrum()

s.plot_phi_pert_spectrum()
s.plot_vrt_spectrum()
s.plot_div_spectrum()

s.plot_physical_fields()

