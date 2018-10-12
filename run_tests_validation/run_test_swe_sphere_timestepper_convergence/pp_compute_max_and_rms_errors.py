#! /usr/bin/env python3

import numpy as np
import sys
import math

from sweet.postprocessing.SphereDataPhysicalDiff import *

s = SphereDataPhysicalDiff(sys.argv[1], sys.argv[2])
print(sys.argv[3]+"\t"+str(s.norm_l1_value)+"\t"+str(s.norm_l2_value)+"\t"+str(s.norm_linf_value))

