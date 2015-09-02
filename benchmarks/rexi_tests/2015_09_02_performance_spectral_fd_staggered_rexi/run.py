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

os.chdir('../../../')


if socket.gethostname() == "inwest":
	print "Running on inwest"
	os.environ['OMP_PROC_BIND'] = "TRUE"
	os.environ['OMP_NUM_THREADS'] = "40"
elif socket.gethostname() == "martinium":
	print "Running on martinium"
	os.environ['OMP_PROC_BIND'] = "TRUE"
	os.environ['OMP_NUM_THREADS'] = "4"


subprocess.call('scons --compiler=intel --program=swe_rexi --spectral-space=enable --libfft=enable --spectral-dealiasing=disable --mode=release '.split(' '), shell=False)
subprocess.call('scons --compiler=intel --program=swe_rexi --spectral-space=disable --libfft=enable --spectral-dealiasing=disable --mode=release '.split(' '), shell=False)

binary_spec = './build/swe_rexi_spectral_libfft_intel_release'
if not os.path.isfile(binary_spec):
	print "Binary "+binary_spec+" not found"
	sys.exit(1)

binary_cart = './build/swe_rexi_libfft_intel_release'
if not os.path.isfile(binary_cart):
	print "Binary "+binary_cart+" not found"
	sys.exit(1)

#
# run for 1 seconds
#
max_time = 1

#
# time step size for coarse time steps
#
rexi_dt = 0.1

#
# order of time step for RK
# Use order 4 to make time errors very small to make the spatial error dominate
#
timestep_order = 4

print "Max simulation time: "+str(max_time)
print "Time step size for REXI time step: "+str(rexi_dt)
print "Time step order: "+str(timestep_order)

#
# default params
#
default_params += '-f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -s 1 -t '+str(max_time)

# Use higher-order time stepping?
default_params += ' -R '+str(timestep_order)


cfl=0.05


# FD/Spectral time stepping
run_method_0 = True

# rexi
run_method_rexi = True



# resolutions
#N_list = [16, 32, 64, 128, 256, 512]
N_list = [16, 32, 64, 128, 256]
N_list = [16, 32, 64]

# REXI parameters
h_list = [0.2]
m_list = [256, 512]

def extract_errors(output):
	match_list = [
		'DIAGNOSTICS ANALYTICAL RMS H:',
		'DIAGNOSTICS ANALYTICAL RMS U:',
		'DIAGNOSTICS ANALYTICAL RMS V:',
		'DIAGNOSTICS ANALYTICAL MAXABS H:',
		'DIAGNOSTICS ANALYTICAL MAXABS U:',
		'DIAGNOSTICS ANALYTICAL MAXABS V:',
		'Simulation time (seconds):'
	]

	vals = ["x" for i in range(len(match_list))]

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


header = ['h_rms', 'u_rms', 'v_rms', 'h_max', 'u_max', 'v_max', 'simtime']

###########################################################################
###########################################################################
# FINITE DIFFERENCES
###########################################################################
###########################################################################


#
title = "Finite-differences, spectral space, C-grid, RK"+str(timestep_order)+", CFL="+str(cfl)
#
if run_method_0:
	print
	print "WARNING!!!! THE ERROR IN VELOCITIES IS NOT VALID DUE TO STAGGERED GRID!"
	print
	print "#TI "+title
	print "#TX RMS, L2, convergence"
	print "#TY Resolution NxN"
	print "N	"+"\t".join(header)+"\t"+"\t".join(['conv('+h+')' for h in header])
	sys.stdout.flush()
	prev_vals = [0.0 for i in range(6)]
	for n in N_list:
		command = binary_spec+' '+default_params
		command += ' -C '+str(cfl)
		command += ' -S 0'
		command += ' --timestepping-mode 0'
		command += ' --staggering 1'
		command += ' -N '+str(n)
		# WARNING: Here we add some viscosity!!!
		# compute hyperviscosity
		# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
		command += ' -U '+str(hyperviscosity[n])

		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		if err != '':
			print "******************"
			print err
			print "******************"
		vals = extract_errors(output)

		conv_rate = [str(float(prev_vals[i])/float(vals[i])) if float(vals[i]) != 0 else 'x' for i in range(6)]
		prev_vals = vals[:]

		print str(n)+"\t"+"\t".join(vals)+"\t"+"\t".join(conv_rate)
		sys.stdout.flush()

#
title = "Finite-differences, spectral space, A-grid, RK"+str(timestep_order)+", CFL="+str(cfl)
#
if run_method_0:
	print
	print "#TI "+title
	print "#TX RMS, L2, convergence"
	print "#TY Resolution NxN"
	print "N	"+"\t".join(header)+"\t"+"\t".join(['conv('+h+')' for h in header])
	sys.stdout.flush()
	prev_vals = [0.0 for i in range(6)]
	for n in N_list:
		command = binary_spec+' '+default_params
		command += ' -C '+str(cfl)
		command += ' -S 0'
		command += ' --timestepping-mode 0'
		command += ' --staggering 0'
		command += ' -N '+str(n)
		# WARNING: Here we add some viscosity!!!
		# compute hyperviscosity
		# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
		command += ' -U '+str(hyperviscosity[n])

		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		if err != '':
			print "******************"
			print err
			print "******************"
		vals = extract_errors(output)

		conv_rate = [str(float(prev_vals[i])/float(vals[i])) if float(vals[i]) != 0 else 'x' for i in range(6)]
		prev_vals = vals[:]

		print str(n)+"\t"+"\t".join(vals)+"\t"+"\t".join(conv_rate)
		sys.stdout.flush()

#
title = "Finite-differences, Cartesian space, C-grid, RK"+str(timestep_order)+", CFL="+str(cfl)
#
if run_method_0:
	print
	print "WARNING!!!! THE ERROR IN VELOCITIES IS NOT VALID DUE TO STAGGERED GRID!"
	print
	print "#TI "+title
	print "#TX RMS, L2, convergence"
	print "#TY Resolution NxN"
	print "N	"+"\t".join(header)+"\t"+"\t".join(['conv('+h+')' for h in header])
	sys.stdout.flush()
	prev_vals = [0.0 for i in range(6)]
	for n in N_list:
		command = binary_cart+' '+default_params
		command += ' -C '+str(cfl)
		command += ' -S 0'
		command += ' --timestepping-mode 0'
		command += ' --staggering 1'
		command += ' -N '+str(n)
		# WARNING: Here we add some viscosity!!!
		# compute hyperviscosity
		# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
		command += ' -U '+str(hyperviscosity[n])

		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		if err != '':
			print "******************"
			print err
			print "******************"

		vals = extract_errors(output)

		conv_rate = [str(float(prev_vals[i])/float(vals[i])) if float(vals[i]) != 0 else 'x' for i in range(6)]
		prev_vals = vals[:]

		print str(n)+"\t"+"\t".join(vals)+"\t"+"\t".join(conv_rate)
		sys.stdout.flush()

#
title = "Finite-differences, Cartesian space, A-grid, RK"+str(timestep_order)+", CFL="+str(cfl)
#
if run_method_0:
	print
	print "#TI "+title
	print "#TX RMS, L2, convergence"
	print "#TY Resolution NxN"
	print "N	"+"\t".join(header)+"\t"+"\t".join(['conv('+h+')' for h in header])
	sys.stdout.flush()
	prev_vals = [0.0 for i in range(6)]
	for n in N_list:
		command = binary_cart+' '+default_params
		command += ' -C '+str(cfl)
		command += ' -S 0'
		command += ' --timestepping-mode 0'
		command += ' --staggering 0'
		command += ' -N '+str(n)
		# WARNING: Here we add some viscosity!!!
		# compute hyperviscosity
		# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
		command += ' -U '+str(hyperviscosity[n])

		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		if err != '':
			print "******************"
			print err
			print "******************"

		vals = extract_errors(output)

		conv_rate = [str(float(prev_vals[i])/float(vals[i])) if float(vals[i]) != 0 else 'x' for i in range(6)]
		prev_vals = vals[:]

		print str(n)+"\t"+"\t".join(vals)+"\t"+"\t".join(conv_rate)
		sys.stdout.flush()



###########################################################################
###########################################################################
# SPECTRAL METHOD SOLVER
###########################################################################
###########################################################################

#
title = "Spectral solver, Spectral space, A-grid, RK"+str(timestep_order)+", CFL="+str(cfl)
#
# spectral time stepping with spectral derivatives
#
if run_method_0:
	print
	print "#TI "+title
	print "#TX RMS, L2, convergence"
	print "#TY Resolution NxN"
	print "N	"+"\t".join(header)+"\t"+"\t".join(['conv('+h+')' for h in header])
	sys.stdout.flush()
	prev_vals = [0.0 for i in range(6)]
	for n in N_list:
		command = binary_spec+' '+default_params
		command += ' -C '+str(cfl)
		command += ' -S 1'		# activate derivatives in spectral space
		command += ' --timestepping-mode 0'
		command += ' --staggering 0'
		command += ' -N '+str(n)
		# WARNING: Here we add some viscosity!!!
		# compute hyperviscosity
		# http://math.boisestate.edu/~wright/research/FlyerEtAl2012.pdf
		command += ' -U '+str(hyperviscosity[n])

		#print command
		p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
		output, err = p.communicate()

		if err != '':
			print "******************"
			print err
			print "******************"

		vals = extract_errors(output)

		conv_rate = [str(float(prev_vals[i])/float(vals[i])) if float(vals[i]) != 0 else 'x' for i in range(6)]
		prev_vals = vals[:]

		print str(n)+"\t"+"\t".join(vals)+"\t"+"\t".join(conv_rate)
		sys.stdout.flush()


###########################################################################
###########################################################################
# REXI
###########################################################################
###########################################################################

#
#
# spectral time stepping with spectral derivatives
#
if run_method_rexi:
	for rexi_h in h_list:
		for rexi_m in m_list:
			title = "REXI solver (h="+str(rexi_h)+", m="+str(rexi_m)+"), A-grid, DT="+str(-rexi_dt)
			print
			print "#TI "+title
			print "#TX RMS, L2, convergence"
			print "#TY Resolution NxN"
			print "N	"+"\t".join(header)+"\t"+"\t".join(['conv('+h+')' for h in header])
			sys.stdout.flush()
			prev_vals = [0.0 for i in range(6)]
			for n in N_list:
				command = binary_cart+' '+default_params
				command += ' -C '+str(-rexi_dt)
				command += ' --timestepping-mode 1'	# REXI
				command += ' -S 0'		# deactivate derivatives in spectral space for DataArray
				command += ' --staggering 0'
				command += ' --rexi-h '+str(rexi_h)
				command += ' --rexi-m '+str(rexi_m)
				command += ' -N '+str(n)

				#print command
				p = subprocess.Popen(command.split(' '), stdout=PIPE, stderr=PIPE, env=os.environ)
				output, err = p.communicate()

				if err != '':
					print "******************"
					print err
					print "******************"

				vals = extract_errors(output)

				conv_rate = [str(float(prev_vals[i])/float(vals[i])) if float(vals[i]) != 0 else 'x' for i in range(6)]
				prev_vals = vals[:]

				print str(n)+"\t"+"\t".join(vals)+"\t"+"\t".join(conv_rate)
				sys.stdout.flush()



print("FIN")
