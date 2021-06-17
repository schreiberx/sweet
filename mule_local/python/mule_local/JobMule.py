import os
import sys

import matplotlib
d = os.getenv('DISPLAY')
if d == None or not d in [':0', ':0.0', ':1']:
    matplotlib.use('agg')

from mule.InfoError import *

from mule_local.JobGeneration import *
from mule_local.JobCompileOptions import *
from mule_local.JobRuntimeOptions import *

from mule.JobPlatforms import *
from mule.JobPlatformResources import *
from mule.JobParallelization import *
from mule.JobParallelizationDimOptions import *


__all__ = ['JobPlatforms', 'JobPlatformResources', 'InfoError', 'JobGeneration', 'JobCompileOptions', 'JobRuntimeOptions', 'JobParallelization', 'JobParallelizationDimOptions']


if __name__ == "__main__":
    if os.getenv("MULE_SOFTWARE_ROOT") is None:
        print("*"*80)
        print("* MULE_SOFTWARE_ROOT environment variable not set")
        print("*")
        print("* Please load SWEET environment first!")
        print("*")
        print("* See README in root folder")
        print("*"*80)
        sys.exit(1)

    p = SWEETJobGeneration()
