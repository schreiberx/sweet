#! /usr/bin/env python3

import sys
import math
import glob

from SWEET import *
from SWEETPostprocessingJobsData import *
from SphereDataPhysicalDiff import *

#from sweet.postprocessing.SphereDataPhysicalDiff import *

from sweet.postprocessing.pickle_SphereDataPhysicalDiff import *


if __name__ == '__main__':
	p = pickle_SphereDataPhysicalDiff()
	p.run()
