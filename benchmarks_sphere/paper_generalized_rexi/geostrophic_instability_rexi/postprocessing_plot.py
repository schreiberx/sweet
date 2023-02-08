#! /usr/bin/env python3

#from mule import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.plotting.Plotting import *

#sys.path.append('../')
#import pretty_plotting as pp
#sys.path.pop()


#
# Outputfile specified?
#
if len(sys.argv) >= 2:
	outfile = sys.argv[1]
else:
	outfile = "output_rexi_time_integration.pdf"



#
# Load data
#
if len(sys.argv) >= 3:
	jd = JobsData(job_dirs=sys.argv[2:])
else:
	jd = JobsData()



#
# Create groups
#
groups = ['runtime.rexi_method', 'runtime.rexi_files_coefficients.0.unique_id_string']
c = JobsDataConsolidate(jd)
job_groups = c.create_groups(groups)

print("Groups:")
for key, g in job_groups.items():
	print(key)


tagname_x = 'runtime.user_defined_parameters.lambda_imag.value'
tagname_y = 'output.error'


def data_filter(x, y, jobdata):
	if float(y) > 1e1:
		return True
	return False


#
# Make ready for plotting
#
d = JobsData_GroupsPlottingScattered(
		job_groups,
		tagname_x,
		tagname_y,
		data_filter = data_filter
	)

data_plotting = d.get_data_float()

# Make pretty
for key, data in data_plotting.items():
	#data['label'] = pp.get_pretty_name(key)
	data['label'] = key.replace('_', '\\_')

#
# Plot!
#
p = Plotting_ScatteredData()

def lambda_fun(p):
	p.ax.set_yscale("log", nonposy='clip')
	p.ax.set_ylim(1e-15, 10)
	p.ax.set_xscale("log", nonposx='clip')

	params = {'legend.fontsize': 8,
			  'legend.handlelength': 2}
	plt.rcParams.update(params)
	

p.plot(
	data_plotting = data_plotting,
	xlabel = "Stiffness ($\lambda \Delta t$)",
	ylabel = "Error",
	title = "Single time step, REXI time integration",
	outfile = outfile,
	lambda_fun = lambda_fun,
)

