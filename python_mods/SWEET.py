import os
import sys

import matplotlib
d = os.getenv('DISPLAY')
if d == None or not d in [':0', ':0.0', ':1']:
	matplotlib.use('agg')

from InfoError import *

from SWEETPlatforms import *
from SWEETPlatformResources import *

from SWEETJobGeneration import *
from SWEETCompileOptions import *
from SWEETRuntimeOptions import *

from SWEETParallelization import *
from SWEETParallelizationDimOptions import *


__all__ = ['SWEETPlatforms', 'SWEETPlatformResources', 'InfoError', 'SWEETJobGeneration', 'SWEETCompileOptions', 'SWEETRuntimeOptions', 'SWEETParallelization', 'SWEETParallelizationDimOptions']


if __name__ == "__main__":
	if os.getenv("SWEET_ROOT") is None:
		print("*"*80)
		print("* SWEET_ROOT environment variable not set")
		print("*")
		print("* Please load SWEET environment first!")
		print("*")
		print("* See README in root folder")
		print("*"*80)
		sys.exit(1)

	p = SWEETJobGeneration()
