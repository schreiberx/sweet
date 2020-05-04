#! /usr/bin/env python3

import sys
import math

from mule_local.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D




groups = [
    'runtime.timestepping_method',
    'runtime.timestepping_order',
]


tagnames_y = [
    'output.error_end_linf_h_pert',
]


j = JobsData(verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(" + "+key)



for tagname_y in tagnames_y:
    print("*"*80)
    print("Processing tagname "+tagname_y)
    print("*"*80)

    tagname_x = 'runtime.space_res_physical'

    if True:
    	"""
    	Plotting format
    	"""
    	d = JobsData_GroupsPlottingScattered(
    			job_groups,
    			tagname_x,
    			tagname_y,
    			meta_attribute_name = 'runtime.timestepping_order',
#    			data_filter = data_filter
    		)

    	for group_name, group_data in d.get_data_float().items():
    		print("*"*80)
    		print("Group: "+group_name)

    		if len(group_data['meta_values']) == 0:
    			raise Exception("No data in group found")

    		prev_value = -1.0
    		conv = '-'
    		convergence_order = None
    		for (x, y, convergence_order_) in zip(group_data['x_values'], group_data['y_values'], group_data['meta_values']):

    			if prev_value > 0:
    				#conv = y/prev_value
    				conv = prev_value/y
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


    		# Cubic for SL with cubic spatial interpolation
    		if '_sl_' in group_name:
    			conv_test_range_end = len(group_data['x_values'])
    			conv_test_range_start = conv_test_range_end-4
    			target_conv = 16
    			error_tolerance_convergence = 0.2
    		else:
    			# Convergence order is stored in meta value
    			conv_test_range_end = len(group_data['x_values'])
    			conv_test_range_start = conv_test_range_end-4
    			target_conv = pow(2.0, meta)
    			error_tolerance_convergence = 0.25


    		print(" + range start/end: "+str(conv_test_range_start)+", "+str(conv_test_range_end))
    		print(" + error_tolerance_convergence: "+str(error_tolerance_convergence))

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
    				#conv = y/prev_value
    				conv = prev_value/y
    			elif prev_value == 0:
    				conv = '[error=0]'

    			error_convergence = '-'
    			if isinstance(conv, float):
    				error_convergence = abs(conv - target_conv)/target_conv

    			print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv)+"\terror: "+str(error_convergence))

    			if error_convergence != '-':
    				if error_convergence > error_tolerance_convergence:
    					print("Error: "+str(error_convergence))
    					raise Exception("Convergence exceeds tolerance of "+str(error_tolerance_convergence))

    			prev_value = y

    		print("[OK]")


    	print("*"*80)
    	print("Convergence tests successful")
    	print("*"*80)
