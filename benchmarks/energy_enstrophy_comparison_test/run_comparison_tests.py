#! /usr/bin/python2

import subprocess
import sys
import os
import time

default_params = ''


# 0: radial dam break
# 1: gaussian
# 2: balanced steady state u
# 3: balanced steady state v
# 4: diamond initial condition
default_params += ' -s 1'

N = 512

# domain size
default_params += ' -n '+str(N)+' -m '+str(N)

# coriolis
f_term = 0.0001

# CFL
default_params += ' -C 0.01'

# verbosity
default_params += ' -v 3'

# how many days? (last value)
default_params += ' -t '+str(60*60*24*30)
#default_params += ' -t '+str(10)

# compiler to use
compiler='gnu'

# KMP affinities for long-term runs
os.environ["KMP_AFFINITY"] = "granularity=core,compact"

# OMP affinity
os.environ["OMP_PLACES"] = "cores"



execs = [
	# swe
	['swe',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=disable',
		'./build/example_swe_'+compiler+'_release',
		''
		],

	# spectral differences active
	['swe_spectral',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=enable',
		'./build/example_swe_spectral_'+compiler+'_release',
		'-S 1'
		],

	# staggered grid
	['swe_staggerred',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe_staggered --spectral-space=disable',
		'./build/example_swe_staggered_'+compiler+'_release',
		''
		],

	# swe
	['swe_f',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=disable',
		'./build/example_swe_'+compiler+'_release',
		'-f '+str(f_term)
		],

	# spectral differences active
	['swe_spectral_f',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=enable',
		'./build/example_swe_spectral_'+compiler+'_release',
		'-S 1 -f '+str(f_term)
		],

	# staggered grid
	['swe_staggerred_f',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe_staggered --spectral-space=disable',
		'./build/example_swe_staggered_'+compiler+'_release',
		'-f '+str(f_term)
		],

	# swe
	['swe_lf',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=disable',
		'./build/example_swe_'+compiler+'_release',
		'-F 1'
		],

	# spectral differences active
	['swe_spectral_lf',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=enable',
		'./build/example_swe_spectral_'+compiler+'_release',
		'-S 1 -F 1'
		],

	# staggered grid
	['swe_staggerred_lf',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe_staggered --spectral-space=disable',
		'./build/example_swe_staggered_'+compiler+'_release',
		'-F 1'
		],


	# swe
	['swe_ud',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=disable',
		'./build/example_swe_'+compiler+'_release',
		'-W 1'
		],

	# spectral differences active
	['swe_spectral_ud',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=enable',
		'./build/example_swe_spectral_'+compiler+'_release',
		'-S 1 -W 1'
		],

	# staggered grid
	['swe_staggerred_ud',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe_staggered --spectral-space=disable',
		'./build/example_swe_staggered_'+compiler+'_release',
		'-W 1'
		],




	# swe
	['swe_lf_C0_9',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=disable',
		'./build/example_swe_'+compiler+'_release',
		'-F 1 -C 0.9'
		],

	# spectral differences active
	['swe_spectral_lf_C0_9',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe --spectral-space=enable',
		'./build/example_swe_spectral_'+compiler+'_release',
		'-S 1 -F 1 -C 0.9'
		],

	# staggered grid
	['swe_staggerred_lf_C0_9',
		'scons --compiler='+compiler+' --gui=disable --compile-program=swe_staggered --spectral-space=disable',
		'./build/example_swe_staggered_'+compiler+'_release',
		'-F 1 -C 0.9'
		],
]


curdir_name = os.getcwd()
print ("Current working directory: "+curdir_name)

os.chdir('../../')
#subprocess.call('make clean'.split(' '), shell=False)


for e in execs:
	name=e[0]
	compiler_cmd=e[1]
	binary=e[2]
	params=e[3]+' '+default_params

	command=[binary]
	command.extend(params.split(' '))

	subprocess.call(compiler_cmd.split(' '), shell=False)

	if not os.path.isfile(binary):
		print "Binary "+binary+" not found"
		sys.exit(1)


for e in execs:
	name=e[0]
	compiler_cmd=e[1]
	binary=e[2]
	params=e[3]+' '+default_params

	command=[binary]
	command.extend(params.split(' '))

	print
	print("Executing "+' '.join(command))
	with open(curdir_name+'/output_'+name+'.csv', 'w') as outfile:
		startTime = time.time()
		subprocess.call(command, stdout=outfile, shell=False)
		endTime = time.time()
		outfile.write('Time: '+str(endTime-startTime)+" secs\n")
	print("DONE")
	print


print("FIN")
