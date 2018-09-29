
import matplotlib
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
	p = SWEETJobGeneration()
