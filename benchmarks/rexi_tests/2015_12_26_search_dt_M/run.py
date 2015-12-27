#! /usr/bin/python2

import subprocess
import sys
import os
import time
from subprocess import PIPE
import socket

default_params = ''


output_file_prefix = 'output'
if len(sys.argv) > 1:
	output_file_prefix = sys.argv[1]


#
# http://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python
#
class Logger(object):
	def __init__(self, filename="Default.log"):
		self.terminal = sys.stdout
		self.log = open(filename, "w")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)


#
# SCENARIO
#
# 0: radial dam break
# 1: gaussian
# 2: balanced steady state u
# 3: balanced steady state v
# 4: diamond initial condition
# 5: waves
default_params = ' --initial-freq-x-mul=2 --initial-freq-y-mul=1 '
scenario_name = "Sin*cos waves"


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



subprocess.call('scons --program=swe_rexi --rexi-parallel-sum=enable --spectral-space=enable --spectral-dealiasing=disable --mode=release '.split(' '), shell=False)

binary = './build/swe_rexi_spectral_libfft_rexipar_gnu_release'
if not os.path.isfile(binary):
	print "Binary "+binary+" not found"
	sys.exit(1)

#
# run for 1 seconds
#
#max_time = 1


#
# order of time step for RK
# Use order 4 to make time errors very small to make the spatial error dominate
#
timestep_order = 4

#
# default params
#
default_params += ' -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 '

# Use higher-order time stepping?
default_params += ' -R '+str(timestep_order)




# time step size for coarse time steps
dt_list = [16/pow(2.0, i) for i in range(20, 0, -1)]

# epsilons
#eps_list = [1, 0.1, 0.01]
#eps_list = [0.01, 0.1, 1]
eps_list = [1]

# resolutions
#N_list = [16, 32, 64, 128, 256, 512]
#N_list = [16, 32, 64, 128, 256]
N_list = [32, 64]

# h values for REXI
h_list = [0.2]

# M values for REXI
M_list = []
M = 1
while M < 20000:
	M_list.append(M)
	M *= 2;

# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
hyperviscosity = {}
for n in N_list:
	hyperviscosity[n] = 4.*pow(float(n), float(-4))
	hyperviscosity[n] = 0


print "Time step size for coarse time step: "+str(dt_list)
print "Time step order: "+str(timestep_order)
print "Search range for h: "+str(h_list)
print "Search range for M: "+str(M_list)
print "Used hyperviscosity: "+str(hyperviscosity)



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


print
print "Running with time stepping mode 1 (search) L_rms:"

for h in h_list:
	for n in N_list:
		for eps in eps_list:

			# redirect output
			filename = curdir_name+'/'+output_file_prefix+"_n"+str(n)+"_eps"+str(eps)+"_h"+str(h)+".csv"
			print "Writing output to "+filename
			fd = open(filename, "w")

			fd.write("#TI res="+str(n)+"x"+str(n)+", eps="+str(eps)+", h="+str(h)+", 10 TS, "+scenario_name+" scenario\n")
			fd.write("#TX REXI parameter M\n")
			fd.write("#TY size of time step\n")
			fd.write("dT\M")
			for M in M_list:
				fd.write("\t"+str(M))
			fd.write("\n")
			fd.flush()

			for dt in dt_list:
				fd.write(str(dt))
				fd.flush()

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
					command += ' -t '+str(dt)
					command += ' -S 1'
					command += ' --use-specdiff-for-complex-array 1 '

					p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
					output, err = p.communicate()

					vals = extract_errors(output)
					fd.write("\t"+str(vals[0]))
					fd.flush()

					print command
					print vals
				fd.write("\n")


print("FIN")
