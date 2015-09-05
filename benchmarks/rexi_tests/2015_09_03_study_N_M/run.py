#! /usr/bin/python2

import subprocess
import sys
import os
import time
from subprocess import PIPE
import socket

default_params = ''


# 0: radial dam break
# 1: gaussian
# 2: balanced steady state u
# 3: balanced steady state v
# 4: diamond initial condition
default_params += ' -s 1 '


curdir_name = os.getcwd()
print ("Current working directory: "+curdir_name)

os.chdir('../../../')


if socket.gethostname() == "inwest":
	print "Running on inwest"
	os.environ['OMP_PROC_BIND'] = "TRUE"
	os.environ['OMP_NUM_THREADS'] = "10"
elif socket.gethostname() == "martinium":
	print "Running on martinium"
	os.environ['OMP_PROC_BIND'] = "TRUE"
	os.environ['OMP_NUM_THREADS'] = "4"



subprocess.call('scons --program=swe_rexi --spectral-space=enable --spectral-dealiasing=disable --mode=release '.split(' '), shell=False)

binary = './build/swe_rexi_spectral_libfft_omp_gnu_release'
if not os.path.isfile(binary):
	print "Binary "+binary+" not found"
	sys.exit(1)


#
# run for 1 seconds
#
max_time = 1

#
# time step size for coarse time steps
#
dt = 0.1

#
# order of time step for RK
# Use order 4 to make time errors very small to make the spatial error dominate
#
timestep_order = 4

print "Max simulation time: "+str(max_time)
print "Time step size for coarse time step: "+str(dt)
print "Time step order: "+str(timestep_order)

#
# default params
#
default_params += ' -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -s 1 -t '+str(max_time)

# Use higher-order time stepping?
default_params += ' -R '+str(timestep_order)



# FD/Spectral time stepping
run_method_0 = False

# extensive search
run_method_1 = False
# small search
run_method_1_search = True

# analytical solution vs. analytical solution
run_method_2 = False



# epsilons
#eps_list = [1, 0.1, 0.01]
#eps_list = [0.01, 0.1, 1]
eps_list = []
eps = 1.0/pow(2.0, 9)
while eps <= 2.0:
	eps_list.append(eps)
	eps *= 2.0;

# resolutions
N_list = [16, 32, 64, 128, 256, 512]

# h values for REXI
h_list = []
h = 1.0/pow(2.0, 9)
while h <= 2.0:
	h_list.append(h)
	h *= 2.0;
h_list = [0.2]

# M values for REXI
M_list = []
M = 1
while M < 50000:
	M_list.append(M)
	M *= 2;



if run_method_1 or run_method_1_search:
	print "Search range for h:"
	print h_list

	print "Search range for M:"
	print M_list



def extract_errors(output):
	match_list = [
		'DIAGNOSTICS ANALYTICAL RMS H:',
		'DIAGNOSTICS ANALYTICAL RMS U:',
		'DIAGNOSTICS ANALYTICAL RMS V:',
		'DIAGNOSTICS ANALYTICAL MAXABS H:',
		'DIAGNOSTICS ANALYTICAL MAXABS U:',
		'DIAGNOSTICS ANALYTICAL MAXABS V:'
	]

	vals = ["x" for i in range(6)]

	ol = output.splitlines(True)
	for o in ol:
		o = o.replace('\n', '')
		o = o.replace('\r', '')
		for i in range(0, len(match_list)):
			m = match_list[i]
			if o[0:len(m)] == m:
				vals[i] = o[len(m)+1:]

	return vals



hyperviscosity = {}
# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
for n in N_list:
	hyperviscosity[n] = 4.*pow(float(n), float(-4))
	hyperviscosity[n] = 0

print "Used hyperviscosity:"
print hyperviscosity

#
# TIME STEPPING MODE 1
#
if run_method_1_search:
	print
	print "Running with time stepping mode 1 (search) L_rms:"

	for eps in eps_list:
		for h in h_list:
			print
			print "#TI Creating study with eps="+str(eps)+" and for h="+str(h)
			print "#TX REXI parameter M"
			print "#TY Resolution NxN"
			
			sys.stdout.write("N\M")
			for M in M_list:
				sys.stdout.write("\t"+str(M))
			sys.stdout.write("\n")

			for n in N_list:
				sys.stdout.write(str(n))
				sys.stdout.flush()

				for M in M_list:
					command = binary+' '+default_params
					command += ' -C '+str(-dt)
					command += ' --timestepping-mode 1'
					command += ' -N '+str(n)
					command += ' --rexi-h '+str(h)
					command += ' --rexi-m '+str(M)
					command += ' -g '+str(eps)
					command += ' -H '+str(eps)
					command += ' -f '+str(eps)

					p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
					output, err = p.communicate()

					sh = 'DIAGNOSTICS ANALYTICAL DIFF H'
					su = 'DIAGNOSTICS ANALYTICAL DIFF U'
					sv = 'DIAGNOSTICS ANALYTICAL DIFF V'

					dh = -1;
					du = -1;
					dv = -1;

					vals = extract_errors(output)
					sys.stdout.write("\t"+str(vals[0]))
					sys.stdout.flush()
				sys.stdout.write("\n")



print("FIN")
