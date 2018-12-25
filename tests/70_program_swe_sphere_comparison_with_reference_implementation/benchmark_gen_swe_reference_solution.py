#! /usr/bin/env python3

## "non-linear barotropically unstable shallow water test case"
## example provided by Jeffrey Whitaker
## https://gist.github.com/jswhit/3845307
##
## Running the script should pop up a window with this image:
## http://i.imgur.com/ZlxR1.png


import numpy as np
import shtns
import sys
import stat
import os

#debug = 10
debug = 0


def outInfo(string, var):
	print(string+": "+str(var.min())+", "+str(var.max()))

class Spharmt(object):
	"""
	wrapper class for commonly used spectral transform operations in
	atmospheric models.  Provides an interface to shtns compatible
	with pyspharm (pyspharm.googlecode.com).
	"""
	def __init__(self,nlons,nlats,ntrunc,rsphere,gridtype='gaussian'):
		"""initialize
		nlons:  number of longitudes
		nlats:  number of latitudes"""
		self._shtns = shtns.sht(ntrunc, ntrunc, 1, shtns.sht_orthonormal+shtns.SHT_NO_CS_PHASE)
		if gridtype == 'gaussian':
			self._shtns.set_grid(nlats,nlons,shtns.sht_quick_init|shtns.SHT_PHI_CONTIGUOUS,0)
		elif gridtype == 'regular':
			self._shtns.set_grid(nlats,nlons,shtns.sht_reg_dct|shtns.SHT_PHI_CONTIGUOUS,0)
		self.lats = np.arcsin(self._shtns.cos_theta)
		self.lons = (2.*np.pi/nlons)*np.arange(nlons)
		self.nlons = nlons
		self.nlats = nlats
		self.ntrunc = ntrunc
		self.nlm = self._shtns.nlm
		self.degree = self._shtns.l
		self.lap = -self.degree*(self.degree+1.0).astype(np.complex)
		self.invlap = np.zeros(self.lap.shape, self.lap.dtype)
		self.invlap[1:] = 1./self.lap[1:]
		self.rsphere = rsphere
		self.lap = self.lap/rsphere**2
		self.invlap = self.invlap*rsphere**2
		
		print("N: "+str(self.nlons)+", "+str(self.nlats))
		print("Mtrunc: "+str(self.ntrunc))
		print("Nlm: "+str(self.nlm))

	def grdtospec(self,data):
		"""compute spectral coefficients from gridded data"""
		return self._shtns.analys(data)
	def spectogrd(self,dataspec):
		"""compute gridded data from spectral coefficients"""
		return self._shtns.synth(dataspec)
	def getuv(self,vrtspec,divspec):
		"""compute wind vector from spectral coeffs of vorticity and divergence"""
		return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)
	def getvrtdivspec(self,u,v):
		"""compute spectral coeffs of vorticity and divergence from wind vector"""
		vrtspec, divspec = self._shtns.analys(u, v)
		return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec





if __name__ == "__main__":
	import time

	# non-linear barotropically unstable shallow water test case
	# of Galewsky et al (2004, Tellus, 56A, 429-440).
	# "An initial-value problem for testing numerical models of the global
	# shallow-water equations" DOI: 10.1111/j.1600-0870.2004.00071.x
	# http://www-vortex.mcs.st-and.ac.uk/~rks/reprints/galewsky_etal_tellus_2004.pdf

	nlons = 196
	nlats = 100
	ntrunc = int(nlons/3)  # spectral truncation (for alias-free computations)

	#if True:
	if False:
		# Debugging with 10 TS
		dt = 60 # time step in seconds
		simtime = 60*10
		output_t = 60
	elif True:
		# Debugging with 10 TS
		# This should be used for the validation test
		dt = 10
		simtime = 60*60*24	# 1 day
		output_t = 60*60
	else:
		#simtime = 66400
		simtime = 60*60*24*6
		dt = 60 # time step in seconds

		# output every day
		output_t = 60*60*24*1

	itmax = int(simtime/dt) # integration length

	rsphere = 6.37122e6 # earth radius
	omega = 7.292e-5 # rotation rate
	grav = 9.80616 # gravity
	hbar = 10.e3 # resting depth
	umax = 80. # jet speed
	phi0 = np.pi/7.; phi1 = 0.5*np.pi - phi0; phi2 = 0.25*np.pi
	en = np.exp(-4.0/(phi1-phi0)**2)
	alpha = 1./3.; beta = 1./15.
	hamp = 120. # amplitude of height perturbation to zonal jet
	efold = 3.*3600. # efolding timescale at ntrunc for hyperdiffusion
	ndiss = 0 # order for hyperdiffusion



	def outputMinMaxSum(i_data, i_prefix):
		if debug == 0:
			return
		
		vmin = i_data.min()
		vmax = i_data.max()
		vsum = i_data.sum()
		
		foo = i_data.ravel()
		vsuminc = 0
		i = 0
		for f in foo:
			vsuminc += f*i
			i += 1
		
		if len(i_data.shape) == 1:
			vsum = vsum*nlons

		print (i_prefix+" | outputMinMaxSum:")
		print ("		| min: "+str(vmin))
		print ("		| max: "+str(vmax))
		print ("		| sum: "+str(vsum))
		if len(i_data.shape) != 1:
			print ("		| suminc: "+str(vsuminc))
		else:
			print ("		| suminc: x")


	def outputSpecMinMaxSum(i_data, i_prefix):
		if debug == 0:
			return

		vmin = i_data.min()
		vmax = i_data.max()
		vsum = i_data.sum()
		
		if len(i_data.shape) != 1:
			print("ERROR: SPEC DATA NOT 1D")
			sys.exit(1)
			
		vsuminc = 0
		i = 0
		for f in i_data:
			vsuminc += f*i
			i += 1

		print (i_prefix+" | outputSpecMinMaxSum:")
		print ("		| min: "+str(vmin.real)+" "+str(vmin.imag))
		print ("		| max: "+str(vmax.real)+" "+str(vmax.imag))
		print ("		| sum: "+str(vsum.real)+" "+str(vsum.imag))
		print ("		| suminc: "+str(vsuminc.real)+" "+str(vsuminc.imag))

	# setup up spherical harmonic instance, set lats/lons of grid
	x = Spharmt(nlons,nlats,ntrunc,rsphere,gridtype='gaussian')
	lons,lats = np.meshgrid(x.lons, x.lats)

	f = 2.*omega*np.sin(lats) # coriolis
	outputMinMaxSum(f, "f")

	def setup(add_bump = True):

		# zonal jet.
		vg = np.zeros((nlats,nlons),np.float)
		u1 = (umax/en)*np.exp(1./((x.lats-phi0)*(x.lats-phi1)))
		outputMinMaxSum(u1, "u1")
		ug = np.zeros((nlats),np.float)
		ug = np.where(np.logical_and(x.lats < phi1, x.lats > phi0), u1, ug)
		ug.shape = (nlats,1)
		ug = ug*np.ones((nlats,nlons),dtype=np.float) # broadcast to shape (nlats,nlonss)

		# height perturbation.
		hbump = hamp*np.cos(lats)*np.exp(-(lons/alpha)**2)*np.exp(-(phi2-lats)**2/beta)
		outputMinMaxSum(hbump, "hbump")

		# initial vorticity, divergence in spectral space
		vrtspec, divspec =  x.getvrtdivspec(ug,vg)

		# Truncate velocities!
		# This is different to the original Whitaker version
		ug, vg = x.getuv(vrtspec, divspec)

		# create hyperdiffusion factor
	#	hyperdiff_fact = np.exp((-dt/efold)*(x.lap/x.lap[-1])**(ndiss/2))

		# solve nonlinear balance eqn to get initial zonal geopotential,
		# add localized bump (not balanced).
		vrtg = x.spectogrd(vrtspec)

		outputMinMaxSum(ug, "ug")
		outputMinMaxSum(vg, "vg")
		outputMinMaxSum(vg, "f")
		
		tmpg1 = ug*(vrtg+f); tmpg2 = vg*(vrtg+f)
		outputMinMaxSum(tmpg1, "tmpg1")
		outputMinMaxSum(tmpg2, "tmpg2")
		
		tmpspec1, tmpspec2 = x.getvrtdivspec(tmpg1,tmpg2)
		outputMinMaxSum(x.spectogrd(tmpspec1), "tmpspec1")
		
		tmpspec2 = x.grdtospec(0.5*(ug**2+vg**2))
		outputMinMaxSum(x.spectogrd(tmpspec2), "tmpspec2")

		phispec = x.invlap*tmpspec1 - tmpspec2
		outputMinMaxSum(x.spectogrd(phispec), "phispec")

		phig = grav*hbar
		phig += x.spectogrd(phispec)
		
		phispec = x.grdtospec(phig)

		if add_bump:
			phispec += x.grdtospec(grav*hbump)

		return vrtspec, divspec, phispec


	def timestep(vrtspec, divspec, phispec):
		# get vort,u,v,phi on grid
		vrtg = x.spectogrd(vrtspec)
		ug, vg = x.getuv(vrtspec,divspec)
		phig = x.spectogrd(phispec)
		
		# compute tendencies.
		tmpg1 = ug*(vrtg+f); tmpg2 = vg*(vrtg+f)
		
		ddivdtspec, dvrtdtspec = x.getvrtdivspec(tmpg1,tmpg2)
		
		dvrtdtspec *= -1
		tmpg = x.spectogrd(ddivdtspec)
		
		tmpg1 = ug*phig; tmpg2 = vg*phig
		
		tmpspec, dphidtspec = x.getvrtdivspec(tmpg1,tmpg2)
				
		dphidtspec *= -1
		tmpspec = x.grdtospec(phig+0.5*(ug**2+vg**2))
		ddivdtspec += -x.lap*tmpspec
		
		#print(np.max(np.abs(dphidtspec)))
		return(dvrtdtspec, ddivdtspec, dphidtspec)

	# Test for optimal geostrophic balance
	vrtspec, divspec, phispec = setup(False)
	tendencies = timestep(vrtspec, divspec, phispec)

	l1_dvrtdt = np.max(np.abs(tendencies[0]))
	l1_ddivdt = np.max(np.abs(tendencies[1]))
	l1_dphidt = np.max(np.abs(tendencies[2]))
	print("l1 dvrtdt = "+str(l1_dvrtdt))
	print("l1 ddivdt = "+str(l1_ddivdt))
	print("l1 dphidt = "+str(l1_dphidt))

	if l1_dphidt/hbar > 1e-12:
		raise Exception("Error too high!")

	if l1_dvrtdt > 1e-16:
		raise Exception("Error too high!")

	if l1_ddivdt > 1e-12:
		raise Exception("Error too high!")

	vrtspec, divspec, phispec = setup(True)

	# initialize spectral tendency arrays
	ddivdtspec = np.zeros(vrtspec.shape+(3,), np.complex)
	dvrtdtspec = np.zeros(vrtspec.shape+(3,), np.complex)
	dphidtspec = np.zeros(vrtspec.shape+(3,), np.complex)
	nnew = 0; nnow = 1; nold = 2



	t = 0
	os.makedirs("job_benchref_solution", exist_ok=True)
	output_file_name = "job_benchref_solution/output_{:s}_t{:020.8f}.csv"

	# Create dummy run.sh script in this folder for automized job processing
	with open("job_benchref_solution/run.sh", "w") as rfile:
		rfile.write("#!/bin/bash\necho \"Dummy\"")
		os.fchmod(rfile.fileno(), stat.S_IRWXU | stat.S_IRWXG | stat.S_IROTH | stat.S_IXOTH)

	def savefile(data, name, t):
		d = x.spectogrd(data)
		d = np.flip(d, 0)
		np.savetxt(output_file_name.format(name, t/(60*60)), d, delimiter="\t")

	savefile(phispec, "prog_phi", t)
	savefile(vrtspec, "prog_vort", t)
	savefile(divspec, "prog_div", t)

	# time loop.
	time1 = time.clock()
	print("ITMAX: "+str(itmax));
	for ncycle in range(itmax):
		t = ncycle*dt

		maxval = x.spectogrd(phispec).max()

		if True:
			# RK1 time stepping
			(vrtdt, divdt, phidt) = timestep(vrtspec, divspec, phispec)
		else:
			# RK2 time stepping
			(vrtdt, divdt, phidt) = timestep(vrtspec, divspec, phispec)
			(vrtdt, divdt, phidt) = timestep(vrtspec+0.5*dt*vrtdt, divspec+0.5*dt*divdt, phispec+0.5*dt*phidt)

		vrtspec += dt*vrtdt
		divspec += dt*divdt
		phispec += dt*phidt

		t += dt

		if t % output_t == 0:
			print("FILEOUTPUT TIMESTEP "+str(ncycle)+", t="+str(t/(60*60))+", maxval="+str(maxval))
			savefile(phispec, "prog_phi", t)
			savefile(vrtspec, "prog_vort", t)
			savefile(divspec, "prog_div", t)
		else:
			if debug > 0:
				print("TIMESTEP "+str(ncycle)+", t="+str(t/(60*60))+", maxval="+str(maxval))

	savefile(phispec, "prog_phi", t)
	savefile(vrtspec, "prog_vort", t)
	savefile(divspec, "prog_div", t)

	vrtg = x.spectogrd(vrtspec)
	ug,vg = x.getuv(vrtspec,divspec)
	phig = x.spectogrd(phispec)

	outputMinMaxSum(vrtg, "vrtg")
	outputMinMaxSum(ug, "ug")
	outputMinMaxSum(vg, "vg")
	outputMinMaxSum(phig, "phig")

