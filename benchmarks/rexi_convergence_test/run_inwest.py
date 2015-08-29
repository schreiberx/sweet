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
default_params += ' -s 1'


curdir_name = os.getcwd()
print ("Current working directory: "+curdir_name)

os.chdir('../../')


if socket.gethostname() == "inwest":
	print "Running on inwest"
	os.environ['OMP_PROC_BIND'] = "TRUE"
	os.environ['OMP_NUM_THREADS'] = "40"
elif socket.gethostname() == "martinium":
	print "Running on martinium"
	os.environ['OMP_PROC_BIND'] = "TRUE"
	os.environ['OMP_NUM_THREADS'] = "4"



subprocess.call('scons --program=swe_rexi --spectral-space=enable --spectral-dealiasing=disable --mode=release --simd=enable'.split(' '), shell=False)

binary = './build/swe_rexi_spectral_gnu_release'
if not os.path.isfile(binary):
	print "Binary "+binary+" not found"
	sys.exit(1)

# run for 3 seconds
max_time = 3

# time step size for 
dt = 0.1

# default params
default_params = '-f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -s 1 -t '+str(max_time)

# Use higher-order time stepping?
#default_params += ' -R 4'

#default_params = '1 -f 1  -g 1 -H 1 -N 64 -X 1 -Y 1  --rexi-h 2 --rexi-m 1 --timestepping-mode 1 --compute-error 1 -v 1  -s 1'


run_method_0 = True
run_method_1 = True
run_method_2 = True

run_method_1_search = True
# resolution for search
N_search = 128

# resolutions
#N_list = [16, 32, 64, 128, 256, 512]
N_list = [16, 32, 64, 128, 256]

# h values for REXI
h_list = []
h = 1.0/pow(2.0, 9)
while h <= 2.0:
	h_list.append(h)
	h *= 2.0;

# M values for REXI
M_list = []
M = 1
while M < 2048:
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

	vals = ["x" for i in range(0, 6)]

	ol = output.splitlines(True)
	for o in ol:
		o = o.replace('\n', '')
		o = o.replace('\r', '')
		for i in range(0, len(match_list)):
			m = match_list[i]
			if o[0:len(m)] == m:
				vals[i] = o[len(m)+1:]

	return vals
	


#
# TIME STEPPING MODE 0
#
# spectral time stepping
# convergence test for resolution and
#
if run_method_0:
	print
	print "Running with time stepping mode 0:"
	print "N	h_rms	u_rms	v_rms	h_max	u_max	v_max"
	sys.stdout.flush()
	for n in N_list:
		command = binary+' '+default_params
		command += ' -C 0.1'
		command += ' --timestepping-mode 0'
		command += ' -N '+str(n)
		# WARNING: Here we add some viscosity!!!
		command += ' -P 0.00000001'

		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		vals = extract_errors(output)
		print str(n)+"\t"+"\t".join(vals)
		sys.stdout.flush()



#
# TIME STEPPING MODE 1
#
if run_method_1:
	print
	print "Running with time stepping mode 1:"
	print "N	rh	rM	h_rms	u_rms	v_rms	h_max	u_max	v_max"
	sys.stdout.flush()

	for h in h_list:
		for M in M_list:
			for n in N_list:
				command = binary+' '+default_params
				command += ' -C '+str(-dt)
				command += ' --timestepping-mode 1'
				command += ' -N '+str(n)
				command += ' --rexi-h '+str(h)
				command += ' --rexi-m '+str(M)

				p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
				output, err = p.communicate()

				sh = 'DIAGNOSTICS ANALYTICAL DIFF H'
				su = 'DIAGNOSTICS ANALYTICAL DIFF U'
				sv = 'DIAGNOSTICS ANALYTICAL DIFF V'

				dh = -1;
				du = -1;
				dv = -1;

				vals = extract_errors(output)
				print str(n)+"\t"+str(h)+"\t"+str(M)+"\t"+"\t".join(vals)
				sys.stdout.flush()



#
# TIME STEPPING MODE 1
#
if run_method_1_search:
	print
	print "Running with time stepping mode 1 (search) L_rms:"

	
	sys.stdout.write("h\M")
	for M in M_list:
		sys.stdout.write("\t"+str(M))
	sys.stdout.write("\n")

	n = N_search
	for h in h_list:
		sys.stdout.write(str(h))
		sys.stdout.flush()

		for M in M_list:
			command = binary+' '+default_params
			command += ' -C '+str(-dt)
			command += ' --timestepping-mode 1'
			command += ' -N '+str(n)
			command += ' --rexi-h '+str(h)
			command += ' --rexi-m '+str(M)

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



#
# TIME STEPPING MODE 2
#
# time stepping with analytical solution
#
if run_method_2:
	print
	print "Running with time stepping mode 2:"
	print "N	h_rms	u_rms	v_rms	h_max	u_max	v_max"
	for n in N_list:
		command = binary+' '+default_params
		command += ' -C '+str(-dt)
		command += ' --timestepping-mode 2'
		command += ' -N '+str(n)

		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		sh = 'DIAGNOSTICS ANALYTICAL DIFF H'
		su = 'DIAGNOSTICS ANALYTICAL DIFF U'
		sv = 'DIAGNOSTICS ANALYTICAL DIFF V'

		h = -1;
		u = -1;
		v = -1;

		vals = extract_errors(output)
		print str(n)+"\t"+"\t".join(vals)
		sys.stdout.flush()



print("FIN")
