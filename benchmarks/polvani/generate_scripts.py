#! /usr/bin/python2

import subprocess
import sys
import os
import time
from subprocess import PIPE
import socket

curdir_name = os.getcwd()
print ("Current working directory: "+curdir_name)

# CFL parameter
CFL=0.3


os.chdir('../../')


output_file_prefix = "output_"


#
# run for 1 seconds
#
max_time = 1000


#
# order of time step for RK
# Use order 4 to make time errors very small to make the spatial error dominate
#
timestep_order = 4



cases_list = [
#		case	filename		R	F	f			g
	[	"E",	"E%Nx%N0000.0000",	.25,	.20,	351.85837720205683,	24.999999254941958	],
]


#
# default params
#


# resolutions
N_list = [64, 128, 256]


# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
hyperviscosity = {}
for n in N_list:
	hyperviscosity[n] = 4.*pow(float(n), float(-4))
	hyperviscosity[n] = 0


print "Time step order: "+str(timestep_order)
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

	il = output.splitlines(True)
	for o in ol:
		o = o.replace('\n', '')
		o = o.replace('\r', '')
		for i in range(0, len(match_list)):
			m = match_list[i]
			if o[0:len(m)] == m:
				vals[i] = o[len(m)+1:]

	return vals


for n in N_list:
	print
	print "Creating study with resolution "+str(n)+"x"+str(n)
	for c in cases_list:
		caseid = c[0]
		filename = c[1]
		param_R = float(c[2])
		param_B = float(c[3])

		param_f = float(c[4])
		param_g = float(c[5])
		param_h = 1

		polvani_case_file = "../../tmp/polvani/Case"+c[0]+"/"+filename.replace("%N", str(n))

		polvani_case_files="BINARY;"
		polvani_case_files+=polvani_case_file+".h;"
		polvani_case_files+=polvani_case_file+".u;"
		polvani_case_files+=polvani_case_file+".v"

		# Generate unique polvani id
		polvani_id = output_file_prefix+"n"+str(n)+"_case"+str(caseid)+"_R"+str(param_R)+"_B"+str(param_B)
		output_filename = curdir_name+'/'+polvani_id+".txt"
		script_filename = curdir_name+'/'+polvani_id+".sh"

		print "Writing batch "+script_filename
		fd = open(script_filename, "w")

		binary = './build/swe_nonstaggered_vector_invariant_spectral_libfft_omp_gnu_release'
		command = binary+' '

		# unit domain
		command += ' -X 1 -Y 1 '

		# Use higher-order time stepping?
		command += ' -C '+str(CFL)
		command += ' -R '+str(timestep_order)
		command += ' -N '+str(n)
		command += ' -g '+str(param_g)
		command += ' -H '+str(param_h)
		command += ' -f '+str(param_f)
		command += ' -t '+str(max_time)
		command += ' -i "'+polvani_case_files+'"'
		command += ' -v 2 '

		fd.write("#! /usr/bin/sh\n")
		fd.write("\n")


		if socket.gethostname() == "inwest":
			fd.write("echo \"Running on inwest\"\n")
			fd.write("export OMP_PROC_BIND=TRUE\n")
			fd.write("export OMP_NUM_THREADS=10\n")

		elif socket.gethostname() == "martinium":
			fd.write("echo \"Running on martinium\"\n")
			fd.write("export OMP_PROC_BIND=TRUE\n")
			fd.write("export OMP_NUM_THREADS=4\n")

		fd.write("""
cd ../../

scons --program=swe_nonstaggered_vector_invariant --spectral-space=enable --spectral-dealiasing=disable --threading=omp --mode=release

"""+command+""" > """+output_filename+"""
""")

