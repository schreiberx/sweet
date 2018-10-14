#! /usr/bin/env python3

from SWEET import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.plotting.Plotting import *

sys.path.append('../')
import pretty_plotting as pp
sys.path.pop()


#
# Load data
#
j = JobsData('job_bench_*', verbosity=0)


#
# Create groups
#
groups = ['runtime.timestepping_method']
c = JobsDataConsolidate(j)
job_groups = c.create_groups(groups)

print("Groups:")
for key, g in job_groups.items():
	print(key)


tagname_x = 'parallelization.num_threads_per_rank'
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

# Make pretty
for key, data in data_plotting.items():
	data['label'] = pp.get_pretty_name(key)

#
# Plot!
#
p = Plotting_ScatteredData()
p.plot(
	data_plotting = data_plotting,
	xlabel = "Number of threads",
	ylabel = "Wallclock time (seconds)",
	title = "Wallclock time",
	outfile = "output_threads_vs_wallclock_time.pdf"
)



#
# Scalability
#
for key, values in data_plotting.items():
	label = key
	x_values = values['x_values']
	y_values = values['y_values']

	# Basis for scalability (number of cores)
	basis_scalability = 1.0

	# Get index of x value for scalability
	i = x_values.index(basis_scalability)
	if i == None:
		raise Exception("Scalability basis not found")

	# Convert to scalability
	values['y_values'] = [y_values[i]/y for y in y_values]

p.plot(
	data_plotting,
	xlabel="Number of threads",
	ylabel="Scalability",
	title = "Scalability",
	outfile="output_threads_vs_scalability.pdf"
)

