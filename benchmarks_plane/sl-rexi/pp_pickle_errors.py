#! /usr/bin/env python3

import sys

from mule.postprocessing.pickle_PlaneDataPhysicalDiff import *
from mule.utils import exec_program


if len(sys.argv) > 1:
	pickle_PlaneDataPhysicalDiff(
				job_dirs=sys.argv[1:],
				#params=["interpolate_ref_to_cmp", "ignore_missing_file"]
				params=["spectral_ref_to_cmp", "ignore_missing_file"]
		)

else:
	pickle_PlaneDataPhysicalDiff(
				#params=["interpolate_ref_to_cmp", "ignore_missing_file"]
				params=["spectral_ref_to_cmp", "ignore_missing_file"]
		)
