
import re
import sys

def matchprefix(l, i_prefix : str):
	return l[:len(i_prefix)] == i_prefix

def getdata(l, i_prefix):
	"""
	return data after prefix or None if prefix was not found
	"""
	if not matchprefix(l, i_prefix):
		return None

	return l[len(i_prefix):]

def postprocessing_read_output_benchmark(input_filename):
	"""
	Read benchmark data from input_filename
	"""
	with open(input_filename, 'r') as f:
		lines = f.readlines()

	solver_groups = {}
	solver_group = None

	for l in lines:
		if l[-1] == '\n':
			l = l[:-1]

		# Plot and restart
		d = getdata(l, 'DATA: solver_group.name: ')
		if d != None:
			solver_groups[d] = {}
			solver_group = solver_groups[d]
			solver_group['name'] = d
			solver_group['simdata'] = []
			continue

		d = getdata(l, 'DATA: solver_group.conv_order: ')
		if d != None:
			solver_group['conv_order'] = d
			continue
			
		d = getdata(l, 'DATA: solver_group.ref_solution: ')
		if d != None:
			solver_group['ref_solution'] = d
			continue

		d = getdata(l, 'INFO: ')
		if d != None:
			# Skip info fields
			continue

		d = getdata(l, 'DATA: ')
		if d != None:
			if not matchprefix(d, 'script_'):
				raise Exception("Expected 'script_' for line\n"+l)

			d = d.split("\t")
			if len(d) != 5:
				raise Exception("ERROR: should have 5 columns")

			# skip invalid nan's
			if d[1] == 'nan':
				print("NaN value detected (probably instability), skipping result")
				continue

			# Start with empty dictionary
			simdata = {}

			ts_name = d[0]

			#
			# extract string identifier of timestepping method for pretty print
			#
			# Example
			# script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_rexi_lc_n_erk_ver1_tso2_tsob2_C000720_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space007_time128
			#
			# script_ was already removed
			#

			# remove script_ prefix
			ts_name = ts_name.replace('script_', '')

			# simparams
			ts_name = re.sub(r".*tsm_", "", ts_name)

			# timestepping orders
			ts_name = re.sub(r"tso._tsob._", "", ts_name)

			# timestep size
			ts_name = re.sub(r"C[0-9]*_", "", ts_name)

			# modes
			ts_name = re.sub(r"M[0-9]{4}_", "", ts_name)

			# REXI stuff
			ts_name = re.sub(r"REXI.*ext[0-9]{2}_", "", ts_name)

			# parallelization
			ts_name = re.sub(r"MPI_.*time[0-9]{3}", "", ts_name)

			if ts_name[-1] == '_':
				ts_name = ts_name[:-1]

			# Setup some data
			simdata['name'] = d[0]
			simdata['ts_name'] = ts_name
			simdata['l1'] = float(d[1])
			simdata['l2'] = float(d[2])
			simdata['linf'] = float(d[3])
			simdata['wallclocktime'] = float(d[4])

			m = re.search('_C([0-9]*)', d[0])
			simdata['dt'] = float(m.group(1))

			solver_group['simdata'].append(simdata)

	return solver_groups
