#! /usr/bin/env python3

import sys
import numpy as np

from mule.postprocessing.JobData import JobData
from mule_local.postprocessing.SphereDataSpectral import SphereDataSpectral
import mule_local.postprocessing.shtnsfiledata as shtnsfiledata
from postprocessing_swe import postprocessing_swe

debug_active = False

print(" ".join(sys.argv))

if len(sys.argv) < 5:
    print("Usage: "+sys.argv[0]+" 'jobdirectory' 'phidatafile.sweet' 'vrtdatafile.sweet' 'divdatafile.sweet'")
    raise Exception("Not enough arguments")

s = postprocessing_swe()
s.setup(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

s.plot_kinetic_energy_distribution()

s.plot_physical_fields()
