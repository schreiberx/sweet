#! /usr/bin/env python3

import os
import sys
import stat
import math

sys.path.append("../../scripts")
from sweet_swe_rexi_plane_and_sphere_params import *
p = sweet_swe_rexi_plane_and_sphere_params()


p.verbosity = 3

p.plane_or_sphere = 'plane'

p.mode_res = -1
p.phys_res = 512

p.simtime = 0.0001
timestep_size_reference = 0.0001
p.timestep_size = timestep_size_reference
p.output_timestep_size = p.timestep_size

p.bench_id = 14

p.rexi_sphere_preallocation = 0

p.g = 1
p.f = 1
p.h = 1
p.h=100
p.domain_size = 1

#p.viscosity = 0.0005
p.viscosity = 0.0

phys_res_list = [16*(2**i) for i in range(0, 7)]

p.nonlinear = 1

# Groups to execute, see below
# l: linear
# ln: linear and nonlinear
#groups = ['l1', 'l2', 'ln1', 'ln2']
groups = ['ln2space']

#if len(sys.argv) < 5:
#	print("Usage: "+str(sys.argv[0])+" [group=l1/l2/ln1/ln2] [tsmethod] [order1] [order2]")
#	sys.exit(1)


if len(sys.argv) > 1:
	groups = [sys.argv[1]]

print("Groups: "+str(groups))

for group in groups:
	# 1st order linear
	if group == 'l1':
		ts_methods = [
			['l_direct',	0,	0,	0,	{'timestep_size': p.simtime}],	# reference solution
			['l_erk',	1,	0],
			['l_irk',	1,	0],
			['l_rexi',	0,	0],
		]

	# 2nd order linear
	if group == 'l2':
		ts_methods = [
			['l_direct',	0,	0,	{'timestep_size': p.simtime}],	# reference solution
			['l_erk',	2,	0],
			['l_cn',	2,	0],
			['l_rexi',	0,	0],
		]

	# 1st order nonlinear
	if group == 'ln1':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['l_erk_n_erk',		1,	1],
			['l_irk_n_erk',		1,	1],
			['ln_erk',		1,	1],
			['l_rexi_n_erk',	1,	1],
		]

	# 1st order nonlinear
	if group == 'ln1test':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['l_erk_n_erk',		1,	1],
			['l_irk_n_erk',		1,	1],
			['ln_erk',		1,	1],
		]

	# 2nd order nonlinear
	if group == 'ln2':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution
			['l_cn_n_erk',		2,	2],
			['l_erk_n_erk',		2,	2],
			['l_irk_n_erk',		2,	2],
			['ln_erk',		2,	2],
			['l_rexi_n_erk',	2,	2],
		]

	# 2nd order nonlinear non-fully-spectral
	if group == 'ln2space':
		ts_methods = [
			['ln_erk',		4,	4],	# reference solution - spectral (128 grid points)
			['ln_erk',		2,	2],	# FD- C-grid
			['l_cn_na_sl_nd_settls', 2,	2],	# SI-SL-SP
#			['l_erk_n_erk',		2,	2],
#			['ln_erk',		2,	2],
#			['l_rexi_n_erk',	2,	2],
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
	if False:
		print("Reference")
		tsm = ts_methods[0]

		p.prefix_string = prefix_string_template+'_ref'
		p.timestepping_method = tsm[0]
		p.timestepping_order = tsm[1]
		p.timestepping_order2 = tsm[2]
		p.phys_res = 512

		if len(tsm) > 4:
			s = tsm[4]
			if 'timestep_size' in s:
				p.timestep_size = s['timestep_size']
		else:
			p.timestep_size = timestep_size_reference

		p.gen_script('script'+p.create_job_id(), 'run.sh')


	for tsm in ts_methods[1:]:

		if group == 'ln2space' and 'ln_erk' in tsm[0]:
			p.staggering = 1
			p.spectralderiv = 0

		if group == 'ln2space' and 'l_cn_na_sl_nd_settls' in tsm[0]:
			p.staggering = 0
			p.spectralderiv = 1

		for phys_res in phys_res_list:

			p.prefix_string = prefix_string_template

			p.timestepping_method = tsm[0]
			p.timestepping_order = tsm[1]
			p.timestepping_order2 = tsm[2]
			p.timestep_size = timestep_size_reference

			p.phys_res = phys_res


			if len(tsm) > 4:
				s = tsm[4]
				p.load_rexi_from_dict(tsm[4])
				if 'timestep_size' in s:
					p.timestep_size = s['timestep_size']

			p.gen_script('script'+p.create_job_id(), 'run.sh')

