#! /usr/bin/env python3

from mule_local.postprocessing.pickle_PlaneDataPhysicalDiff import *
from mule.exec_program import *

pickle_PlaneDataPhysicalDiff(params=["interpolate", "ignore_missing_file"])
