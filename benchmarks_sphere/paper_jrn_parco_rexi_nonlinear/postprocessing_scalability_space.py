#! /usr/bin/env python3

from SWEET import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.plotting.Plotting import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



j = JobsData('job_bench_*', verbosity=100)

c = JobsDataConsolidate(j)
data_table = c.create_data_table_float(
		['runtime.timestepping_method'],
		'parallelization.num_threads_per_rank',
		'output.simulation_benchmark_timings.main_timestepping'
	)

print("Data table:")
c.print_data_table(data_table)
c.write_data_table(data_table, "output_threads_vs_wallclock_time.csv")


data_plotting = c.create_data_plotting_float(
		['runtime.timestepping_method'],
		'parallelization.num_threads_per_rank',
		'output.simulation_benchmark_timings.main_timestepping'
	)


p = Plotting()

#
# Wallclock time
#
p.plot_scattered(data_plotting, xlabel="Number of threads", ylabel="Wallclock time (seconds)", title = "Wallclock time", outfile="output_threads_vs_wallclock_time.pdf")


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

p.plot_scattered(data_plotting, xlabel="Number of threads", ylabel="Scalability", title = "Scalability", outfile="output_threads_vs_scalability.pdf")
