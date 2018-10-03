import matplotlib
matplotlib.use('agg')

import os
import sys

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
