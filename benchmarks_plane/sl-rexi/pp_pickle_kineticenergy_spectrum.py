#! /usr/bin/env python3

from mule.postprocessing.pickle_PlaneData_KineticEnergy import *

# process all time frames
if len(sys.argv) > 1:
	pickle_PlaneData_KineticEnergy(params=[], job_dirs=sys.argv[1:])
else:
	pickle_PlaneData_KineticEnergy(params=[])
