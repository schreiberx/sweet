#! /usr/bin/env python3

from SWEET import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.plotting.Plotting import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



j = JobsData('job_bench_*', verbosity=100)

c = JobsDataConsolidate(j)

d = JobsData_DataTable(
		j,
		'parallelization.num_ranks',
		[
			'output.simulation_benchmark_timings.main_timestepping',
			'output.simulation_benchmark_timings.rexi',
			'output.simulation_benchmark_timings.rexi_setup',
			'output.simulation_benchmark_timings.rexi_shutdown',
			'output.simulation_benchmark_timings.rexi_timestepping',
			'output.simulation_benchmark_timings.rexi_timestepping_solver',
			'output.simulation_benchmark_timings.rexi_timestepping_broadcast',
			'output.simulation_benchmark_timings.rexi_timestepping_reduce',
			'output.simulation_benchmark_timings.rexi_timestepping_miscprocessing',
		]
	)

print("Data table:")

# make nice looking
for i in range(len(d.data[0])):
	d.data[0][i] = d.data[0][i].replace('output.simulation_benchmark_timings.', '')

d.print()


#p.plot_scattered(data_plotting, xlabel="Number of threads", ylabel="Scalability", title = "Scalability", outfile="output_threads_vs_scalability.pdf")
