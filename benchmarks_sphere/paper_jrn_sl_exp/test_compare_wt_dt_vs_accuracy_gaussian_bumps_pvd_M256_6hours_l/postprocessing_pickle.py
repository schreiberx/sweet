#! /usr/bin/env python3

import sys

from mule_local.postprocessing.pickle_SphereDataPhysicalDiff import *

if len(sys.argv) > 1:
    p = pickle_SphereDataPhysicalDiff(sys.argv[1:])
else:
    p = pickle_SphereDataPhysicalDiff()

