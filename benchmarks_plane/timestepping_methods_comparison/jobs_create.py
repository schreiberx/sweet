#! /usr/bin/env python3

import os
import sys
import stat
import math

from sweet_swe_rexi_plane_params import *
p = sweet_swe_rexi_plane_params()



p.verbosity = 2

p.plane_or_sphere = 'plane'

p.mode_res = -1
p.phys_res = 128

p.bench_id = 1

p.rexi_sphere_preallocation = 0

p.g = 1
p.f = 1
p.h = 1
p.domain_size = 1

#p.viscosity = 0.0005
p.viscosity = 0.0

p.simtime = 0.1
p.output_timestep_size = p.simtime

timestep_size_reference = 0.0001
timestep_sizes = [0.0001*(2.0**i) for i in range(0, 11)]


# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
groups = ['l1', 'l2', 'ln1', 'ln2', 'ln4']
#groups = ['ln2test']

#if len(sys.argv) < 5:
#	print("Usage: "+str(sys.argv[0])+" [group=l1/l2/ln1/ln2] [tsmethod] [order1] [order2] [use rexi direct solution]")
#	sys.exit(1)


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
	# 1st order linear
	if group == 'l1':
		ts_methods = [
			['l_direct',	0,	0,	0],	# reference solution
			['l_erk',	1,	0,	0],
			['l_irk',	1,	0,	0],
			['l_rexi',	0,	0,	0],
		]

	# 2nd order linear
	if group == 'l2':
		ts_methods = [
			['l_direct',	0,	0,	0],	# reference solution
			['l_erk',	2,	0,	0],
			['l_cn',	2,	0,	0],
			['l_rexi',	0,	0,	0],
		]

	#	['lg_rexi_lc_erk_nt_sl_nd_erk',
	#	['l_rexi_ns_sl_nd_erk',

	# 1st order nonlinear
	if group == 'ln1':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			['l_erk_n_erk',		1,	1,	0],
			['l_irk_n_erk',		1,	1,	0],
			['ln_erk',		1,	1,	0],
			['l_rexi_n_erk',	1,	1,	0],
		]

	# 2nd order nonlinear
	if group == 'ln2':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			['l_cn_n_erk',		2,	2,	0],
			['l_erk_n_erk',		2,	2,	0],
			['l_irk_n_erk',		2,	2,	0],
			['ln_erk',		2,	2,	0],
			['l_rexi_n_erk',	2,	2,	0],
			['l_rexi_ns_sl_nd_erk',	2,	2,	0],
			['lg_rexi_lc_erk_nt_sl_nd_erk',	2,	2,	0],
		]


	# 4th order nonlinear
	if group == 'ln4':
		ts_methods = [
			['ln_erk',		4,	4,	0],	# reference solution
			#['ln_etdrk',		4,	4,	1],	# reference solution

			['ln_etdrk',		4,	4,	1],
			['ln_erk',		4,	4,	0],
		]



	#
	# OVERRIDE TS methods
	#
	if len(sys.argv) > 4:
		ts_methods = [ts_methods[0]]+[[sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]]


	#
	# add prefix string to group benchmarks
	#
	prefix_string_template = group


	#
	# Reference solution
	#
	if True:
		tsm = ts_methods[0]

		p.prefix_string = prefix_string_template+'_ref'
		p.timestepping_method = tsm[0]
		p.timestepping_order = tsm[1]
		p.timestepping_order2 = tsm[2]
		p.rexi_use_direct_solution = tsm[3]

		if len(tsm) > 4:
			s = tsm[4]
			if 'timestep_size' in s:
				p.timestep_size = s['timestep_size']
		else:
			p.timestep_size = timestep_size_reference

		p.gen_script('script'+p.create_job_id(), 'run.sh')


	for tsm in ts_methods[1:]:
		for p.timestep_size in timestep_sizes:
			p.prefix_string = prefix_string_template

			p.timestepping_method = tsm[0]
			p.timestepping_order = tsm[1]
			p.timestepping_order2 = tsm[2]
			p.rexi_use_direct_solution = tsm[3]

			if len(tsm) > 4:
				s = tsm[4]
				p.load_rexi_from_dict(tsm[4])
				if 'timestep_size' in s:
					p.timestep_size = s['timestep_size']

			p.gen_script('script'+p.create_job_id(), 'run.sh')

