#! /usr/bin/env python3

import sys
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from postprocessing_read_output_benchmark import *

#
# First, use
#   ./postprocessing.py > postprocessing_output.txt
# to generate the .txt file
#

if len(sys.argv) < 3:
	raise Exception("Please specify at least one reference file and files for comparisons")

# input filename for reference data
input_filename_ref = sys.argv[1]

# input filename for comparison
input_filenames = sys.argv[2:]


print("INFO: "+"*"*70)
print("INFO: Input filename reference: "+input_filename_ref)
for i in range(len(input_filenames)):
	print("INFO: Input filename["+str(i)+"]: "+input_filenames[i])
print("INFO: "+"*"*70)

print("INFO: Loading data from reference '"+input_filename_ref+"'")
solver_groups_ref = postprocessing_read_output_benchmark(input_filename_ref)

print("INFO: "+"*"*70)

for input_filename in input_filenames:
	print("INFO: Loading data from '"+input_filename+"'")
	solver_groups = postprocessing_read_output_benchmark(input_filename)
	print("INFO: Comparing...")

	max_error = 1e-10

	for key, solver_group_ref in solver_groups_ref.items():

		print("INFO: solver_group_ref.name: "+solver_group_ref['name'])
		solver_group = solver_groups[key]

		# setup accumulation buffers
		for i in range(len(solver_group_ref['simdata'])):
			simdata_ref = solver_group_ref['simdata'][i]
			simdata = solver_group['simdata'][i]

			err = abs(simdata_ref['l1'] - simdata['l1'])
			print(simdata['name']+"\t"+str(err))

			if err > max_error:
				raise Exception("Maximum error "+str(max_error)+" exceeded")
