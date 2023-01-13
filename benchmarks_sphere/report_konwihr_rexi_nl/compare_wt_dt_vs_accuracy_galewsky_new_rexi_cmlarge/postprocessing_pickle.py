#! /usr/bin/env python3

import sys
import math
import glob

from mule.postprocessing.pickle_SphereDataSpectralDiff import *
from mule.utils import exec_program

# Ugly hack!
#output, retval = exec_program('ls *benchref*/*prog_h* | sort | tail -n 1 | sed "s/.*prog_h//"')
#if retval != 0:
#	print(output)
#	raise Exception("Something went wrong")

#output = output.replace("\n", '')
#output = output.replace("\r", '')

#p = pickle_SphereDataSpectralDiff(output)
p = pickle_SphereDataSpectralDiff()
