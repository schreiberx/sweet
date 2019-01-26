#! /usr/bin/env python3

from mule_local.postprocessing.pickle_PlaneData_KineticEnergy import *

params=["only_last_file", "ignore_missing_file"]

# process only the last time frames (likely the reference time frame)
if len(sys.argv) > 1:
	pickle_PlaneData_KineticEnergy(params=params, job_dirs=sys.argv[1:])
else:
	pickle_PlaneData_KineticEnergy(params=params)
