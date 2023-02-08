#! /usr/bin/env python3

from mule.postprocessing.pickle_PlaneData_Spectrum import *
from mule.utils import exec_program

# process only last file (likely the reference solution)
pickle_PlaneData_Spectrum(params=["only_last_file", "ignore_missing_file"])

# process all files
#pickle_PlaneData_Spectrum(params=[])
