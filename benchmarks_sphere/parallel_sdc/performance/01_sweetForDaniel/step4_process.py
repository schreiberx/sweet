#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.plotting.Plotting import *

sys.path.append('../')
sys.path.pop()

#
# Load data
#
j = JobsData('job_bench_*', verbosity=0)


#
# Create groups
#
groups = ['parallelization.pType']
c = JobsDataConsolidate(j)
job_groups = c.create_groups(groups)

print("Groups:")
for key, g in job_groups.items():
	print(key)


tagname_x = 'parallelization.nProc'
tagname_y = 'output.simulation_benchmark_timings.main_timestepping'

#
# Make ready for plotting
#
d = JobsData_GroupsPlottingScattered(
		job_groups,
		tagname_x,
		tagname_y
	)

data_plotting = d.get_data_float()
print(data_plotting)


#
# Plot!
#

# Utility function
iSym = [0]
lSym = ['o', '^']
def sym():
	s = lSym[iSym[0]]
	iSym[0] += 1
	return s


tBase = data_plotting['Space Parallel']['y_values'][0]
for group, data in data_plotting.items():

	nProc = np.array(data['x_values'])
	tComp = np.array(data['y_values'])
	s = sym()

	plt.figure('wallclock')
	plt.loglog(nProc, tComp, s+'-', label=group)
	if group == 'Space Parallel':
		plt.loglog(nProc, tBase/nProc, '--', c='gray')

	plt.figure('speedup')
	plt.plot(nProc, tBase/tComp, s+'-', label=group)
	if group == 'Space Parallel':
		plt.plot(nProc, nProc, '--', c='gray')
	plt.ylim(0, 10)
	

for figName in ['wallclock', 'speedup']:
	plt.figure(figName)
	plt.xlabel('nProc')
	plt.ylabel(figName)
	plt.legend()
	plt.grid()
	plt.tight_layout()
	plt.savefig(f'output_{figName}.pdf')
