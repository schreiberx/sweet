#! /usr/bin/env python3

import sys
import math

from SWEET import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def get_test_parameters(method, order):
	"""
	Return:
	-------
		( start range, end range, error tolerance )
	"""
	# Test first 4 fields for convergence
	return 0, 4, 0.95



groups = ['runtime.timestepping_method']


tagnames_y = [
	'plane_data_diff_prog_h_pert.res_norm_l1',
#	'plane_data_diff_prog_h_pert.res_norm_linf',
]


j = JobsData('./job_bench_*', verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
	print(" + "+key)

# Filter out errors beyond this value!
def data_filter(x, y, jobdata):

	if y == None:
		return True

	# Filter out NaNs for wallclock time studies
	# NaNs require significantly more computation time
	if math.isnan(y):
		return True

	return False

	if y > 10.0:
		return True

	return False


for tagname_y in tagnames_y:
	print("*"*80)
	print("Processing tagname "+tagname_y)
	print("*"*80)

	tagname_x = 'runtime.timestep_size'

	if True:
		"""
		Plotting format
		"""

		d = JobsData_GroupsPlottingScattered(
				job_groups,
				tagname_x,
				tagname_y,
				meta_attribute_name = 'runtime.timestepping_order',
#				data_filter = data_filter
			)

		for key, group_data in d.get_data_float().items():
			print("*"*80)
			print("Group: "+key)
			prev_value = -1.0
			conv = '-'
			for (x, y, convergence_order) in zip(group_data['x_values'], group_data['y_values'], group_data['meta_values']):
				if prev_value >= 0:
					conv = y/prev_value

				print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv))
				prev_value = y


			print("")
			print("Testing convergence")
			(conv_test_range_start, conv_test_range_end, error_tolerance) = get_test_parameters(key, convergence_order)
			print(" + range start/end: "+str(conv_test_range_start)+", "+str(conv_test_range_end))
			print(" + error_tolerancce: "+str(error_tolerance))

			if len(group_data['meta_values']) < conv_test_range_end:
				raise Exception("Not enough samples to run convergence test")

			for i in range(len(group_data['meta_values'])):
				if group_data['meta_values'][i] != group_data['meta_values'][0]:
					raise Exception("FATAL: Different convergence orders in same test")

			prev_value = -1.0
			conv = '-'
			for i in range(conv_test_range_start, conv_test_range_end):
				x = group_data['x_values'][i]
				y = group_data['y_values'][i]
				meta = group_data['meta_values'][i]

				if prev_value >= 0:
					conv = y/prev_value

				if conv != '-':
					target_conv = pow(2.0, meta)
					error = abs(conv - target_conv)/target_conv

					if error > error_tolerance:
						print("Error: "+str(error))
						raise Exception("Convergence exceeds tolerance of "+str(error_tolerance))
				
				print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv))
				prev_value = y

			print("[OK]")

		print("*"*80)
		print("Convergence tests successful")
		print("*"*80)
"""
		print("Expected convergence: "+str(conv_order))
		print("Measured convergence: "+str(conv_test))
                
		# test these first convergence tests
		if 'ln4' in group_info[0]:
			# Comparing RK4 with 4-th order methods requires a more relaxed test
			# This works and was empirically determined
			test_range = range(2,5)
			max_error_rate = 0.5
		elif 'ln2' in group_info[0]:
			test_range = range(1,3)
			max_error_rate = 0.5
		else:
			test_range = range(1,3)
			max_error_rate = 0.05

		print(test_range)

		for i in test_range:
			if conv_test[i] == 0.0:
				if 'ln4' not in group_info[0]:
					print("Invalid convergence of 0")
					sys.exit(1)

			conv_rate = 2**conv_order
			a = abs(conv_test[i]-conv_rate)/conv_rate
			print("Testing convergence "+str(conv_test[i])+": "+str(a))

		print("	None: No data available")
"""
