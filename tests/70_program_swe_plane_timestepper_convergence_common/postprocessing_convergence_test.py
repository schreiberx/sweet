#! /usr/bin/env python3

import sys
import math

from SWEET import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D





groups = ['runtime.timestepping_method']


tagnames_y = [
	'plane_data_diff_prog_h_pert.res_norm_linf',
	'plane_data_diff_prog_u.res_norm_linf',
	'plane_data_diff_prog_v.res_norm_linf',
]


j = JobsData(verbosity=0)

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
			convergence_order = None
			for (x, y, convergence_order_) in zip(group_data['x_values'], group_data['y_values'], group_data['meta_values']):

				if prev_value > 0:
					conv = y/prev_value
				elif prev_value == 0:
					conv = '[error=0]'

				print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv))
				prev_value = y

				if convergence_order != None:
					if convergence_order != convergence_order_:
						raise Exception("Convergence order mismatch!!!")

				convergence_order = convergence_order_


			print("")
			print("Testing convergence")

			if 'vort' in tagname_y:
				conv_test_range_start = 0
				conv_test_range_end = 4
				error_tolerance = 0.1
			else:
				conv_test_range_start = 0
				conv_test_range_end = 4
				error_tolerance = 0.05

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

				if prev_value > 0:
					conv = y/prev_value
				elif prev_value == 0:
					conv = '[error=0]'

				error = '-'
				if isinstance(conv, float):
					# Convergence order is stored in meta value
					target_conv = pow(2.0, meta)
					error = abs(conv - target_conv)/target_conv

				print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv)+"\terror: "+str(error))

				if error != '-':
					if error > error_tolerance:
						print("Error: "+str(error))
						raise Exception("Convergence exceeds tolerance of "+str(error_tolerance))

				prev_value = y

			print("[OK]")

		print("*"*80)
		print("Convergence tests successful")
		print("*"*80)
