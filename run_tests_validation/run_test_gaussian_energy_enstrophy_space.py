#! /usr/bin/python2

import subprocess
import sys
import os
import time


curdir_name = os.getcwd()
print("Current working directory: "+curdir_name)

os.chdir('../')


sim_time=100
f=0.0001
benchmark_id=2

exec_stuff=[
		# compile command | execution | parameters | superconvergence?
		['--compile-program=swe --spectral-space=enable --spectral-dealiasing=enable ', './build/example_swe_spectral_dealiasing_gnu_release', '-s '+str(benchmark_id)+' -S 1 -C 0.01', True],
		['--compile-program=swe --spectral-space=enable --spectral-dealiasing=disable ', './build/example_swe_spectral_gnu_release', '-s '+str(benchmark_id)+' -S 1 -C 0.01', True],
		['--compile-program=swe --spectral-space=enable --spectral-dealiasing=disable ', './build/example_swe_spectral_gnu_release', '-s '+str(benchmark_id)+' -S 0 -C 0.01', False],
		['--compile-program=swe --spectral-space=disable --spectral-dealiasing=disable ', './build/example_swe_gnu_release', '-s '+str(benchmark_id)+' -S 0 -C 0.01', False],
		['--compile-program=swe_staggered_covariant --spectral-space=enable --spectral-dealiasing=disable ', './build/example_swe_staggered_covariant_spectral_gnu_release', '-s '+str(benchmark_id)+' -S 0 -C 0.01', False],
		['--compile-program=swe_staggered_covariant --spectral-space=disable --spectral-dealiasing=disable ', './build/example_swe_staggered_covariant_gnu_release', '-s '+str(benchmark_id)+' -S 0 -C 0.01', False],
	]

for e in exec_stuff:

	subprocess.call('make clean'.split(' '), shell=False)

	compile_cmd='scons '+e[0]+' --gui=disable --mode=release'
	subprocess.call(compile_cmd.split(' '), shell=False)

	params=" -f "+str(f)+" -v 2 -t "+str(sim_time)+" "+e[2]

	print("COMMAND: "+e[1])

	for R in [1, 2, 3, 4]:
		prev_diff_value=0
		for N in [16, 32, 64, 128, 256]:
			command=e[1]+" "+params+" -N "+str(N)+" -R "+str(R)
			#print("Executing: "+command)
			p = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE)
			result = p.communicate()[0]

			diff_value = -1
			for r in result.split('\n'):
				matchstr="DIAGNOSTICS ENERGY DIFF:"
				#matchstr="DIAGNOSTICS BENCHMARK DIFF H:"
				if r[0:len(matchstr)]==matchstr:
					diff_value=float(r.split('\t')[1])

			if diff_value == -1:
				print "Difference value not found"
				sys.exit(-1)

			if e[3] and False:
				eps=1e-10
				if diff_value > eps:
					print("diff_value "+str(diff_value)+" exceeded superconvergence threshold "+str(eps))
					sys.exit(1)

				conv_rate='S'

			else:
				if diff_value == 0:
					conv_rate = 'X'
				else:
					conv_rate = prev_diff_value/diff_value

			print str(N)+"\t"+str(R)+"\t"+str(diff_value)+"\t"+str(conv_rate)

			prev_diff_value = diff_value

