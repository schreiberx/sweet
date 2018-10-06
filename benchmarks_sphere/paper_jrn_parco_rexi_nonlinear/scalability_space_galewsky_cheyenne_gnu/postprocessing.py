#! /usr/bin/env python3

from SWEET import *
from SWEETPostprocessing import *
from SWEETPlotting import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



j = SWEETPostprocessing(verbosity=0)

data_table = j.create_data_table_float(
		['runtime.timestepping_method'],
		'parallelization.num_threads_per_rank',
		'output.benchmark_timings.main_simulationloop'
	)

print("Data table:")
j.print_data_table(data_table)


data_plotting = j.create_data_plotting_float(
		['runtime.timestepping_method'],
		'parallelization.num_threads_per_rank',
		'output.benchmark_timings.main_simulationloop'
	)


p = SWEETPlotting()

#
# Wallclock time
#
p.plot_scattered(data_plotting, x_label="Number of threads", y_label="Wallclock time (seconds)", title = "Wallclock time", outfile="output_threads_vs_wallclock_time.pdf")
#p.plot_scattered(data_plotting, x_label="Number of threads", y_label="Wallclock time (seconds)", title = "Wallclock time")


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

p.plot_scattered(data_plotting, x_label="Number of threads", y_label="Scalability", title = "Scalability", outfile="output_threads_vs_scalability.pdf")
#p.plot_scattered(data_plotting, x_label="Number of threads", y_label="Scalability", title = "Scalability")
