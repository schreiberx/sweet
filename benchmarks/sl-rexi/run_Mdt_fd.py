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
#default_params = ' --initial-freq-x-mul=2 --initial-freq-y-mul=1 '
default_params = ' '
scenario_name = "Nonlinear-steady"


curdir_name = os.getcwd()
print ("Current working directory: "+curdir_name)

#os.chdir('../../../')

os.environ['OMP_PROC_BIND'] = "TRUE"
os.environ['OMP_NUM_THREADS'] = "8"

subprocess.call('scons --program=swe_rexi --spectral-space=enable --mode=release --threading=off --gui=enable --parareal=none --rexi-parallel-sum=enable --spectral-dealiasing=disable '.split(' '), shell=False)

#subprocess.call('scons --program=swe_rexi --rexi-parallel-sum=enable --spectral-space=enable --spectral-dealiasing=disable --mode=release '.split(' '), shell=False)

binary = './build/swe_rexi_spectral_libfft_gui_rexipar_gnu_release' 
# './build/swe_rexi_spectral_libfft_rexipar_gnu_release'
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
#default_params += ' -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 '
default_params += '-g 1 -f 1 -X 1 -Y 1 -H 10 --nonlinear 1 --compute-error 1 -G 0 -S 0 -v 2 -s 14 --timestepping-mode 1 --staggering 0'

#'-N 128  -C -0.0001  -t 0.001 '

# Use higher-order time stepping?
#default_params += ' -R '+str(timestep_order)

# Max time
T = 0.1

# time step size for coarse time steps
dt_list = [T/pow(2.0, i) for i in range(8, 20, +1)]

# epsilons
#eps_list = [1, 0.1, 0.01]
#eps_list = [0.01, 0.1, 1]
eps_list = [1]

# resolutions
#N_list = [16, 32, 64, 128, 256, 512]
#N_list = [16, 32, 64, 128, 256]
N_list = [ 512]

# h values for REXI
h_list = [0.2]

# M values for REXI
M_list = []
M = 4
while M < 1000:
	M_list.append(M)
	M *= 2;

# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
#hyperviscosity = {}
#for n in N_list:
#	hyperviscosity[n] = 4.*pow(float(n), float(-4))
#	hyperviscosity[n] = 0


print "Time step size for coarse time step: "+str(dt_list)
print "Time step order: "+str(timestep_order)
print "Search range for h: "+str(h_list)
print "Search range for M: "+str(M_list)
#print "Used hyperviscosity: "+str(hyperviscosity)



def extract_errors(output):
	match_list = [
		'DIAGNOSTICS BENCHMARK DIFF H:',
		'DIAGNOSTICS BENCHMARK DIFF U:',
		'DIAGNOSTICS BENCHMARK DIFF V:'
	]

	vals = ["x" for i in range(3)]

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
print "Running with time stepping mode 1:"

for h in h_list:
	for n in N_list:
                # redirect output
                filename = curdir_name+'/'+output_file_prefix+"_n"+str(n)+".csv"
                filename_dump = curdir_name+'/'+output_file_prefix+"_n"+str(n)+"dump.txt"
                print "Writing output to "+filename
                fd = open(filename, "w")
                fd_dump = open(filename_dump, "w")

                fd.write("#TI res="+str(n)+"x"+str(n)+", h="+str(h)+", Max time="+str(T)+" , "+scenario_name+" scenario\n")
                fd.write("#TX REXI parameter M\n")
                fd.write("#TY size of time step\n")
                fd.write("# "+default_params+" \n")
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
                                command += ' -N '+str(n)
                                command += ' --rexi-m '+str(M)
                                command += ' -t '+str(T)

                                p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
                                output, err = p.communicate()
                                fd_dump.write(output)
          
                                vals = extract_errors(output)
                                fd.write("\t"+str(vals[0]))
                                fd.flush()

                                print command
                                print vals
                        fd.write("\n")


print("FIN")
