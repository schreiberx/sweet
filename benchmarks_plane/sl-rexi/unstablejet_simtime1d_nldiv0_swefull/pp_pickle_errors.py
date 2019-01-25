#! /usr/bin/env python3

from mule_local.postprocessing.pickle_PlaneDataPhysicalDiff import *
from mule.exec_program import *

pickle_PlaneDataPhysicalDiff(params=["interpolate_ref_to_cmp", "ignore_missing_file"])
