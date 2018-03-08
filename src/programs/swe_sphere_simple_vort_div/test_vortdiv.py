#! /usr/bin/env python3



import numpy as np
import shtns
import sys


class Spharmt(object):
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
		print("Nlm: "+str(self.nlm))

	def getuv(self,vrtspec,divspec):
		"""compute wind vector from spectral coeffs of vorticity and divergence"""
		return self._shtns.synth((self.invlap/self.rsphere)*vrtspec, (self.invlap/self.rsphere)*divspec)

	def getvrtdivspec(self,u,v):
		"""compute spectral coeffs of vorticity and divergence from wind vector"""
		vrtspec, divspec = self._shtns.analys(u, v)
		return self.lap*self.rsphere*vrtspec, self.lap*rsphere*divspec




if __name__ == "__main__":
	#T = 127
	T = 31
	if True:
	#if False:
		# No dealiasing to avoid truncation of modes
		ntrunc = T
		nlons = (T+1)*2
		nlats = (T+1)
	else:
		nlons = (T+1)*3//2
		nlats = (T+1)*3//2
		ntrunc = T

	rsphere = 6.37122e6 # earth radius
	omega = 7.292e-5 # rotation rate
	grav = 9.80616 # gravity
	hbar = 10.e3 # resting depth
	umax = 80. # jet speed
	phi0 = np.pi/7.; phi1 = 0.5*np.pi - phi0; phi2 = 0.25*np.pi
	en = np.exp(-4.0/(phi1-phi0)**2)
	hamp = 120. # amplitude of height perturbation to zonal jet



	def outputMinMaxSum(i_data, i_prefix):
		vmin = i_data.min()
		vmax = i_data.max()
		
		print (i_prefix+" | outputMinMaxSum:")
		print ("		| min: "+str(vmin))
		print ("		| max: "+str(vmax))

	# setup up spherical harmonic instance, set lats/lons of grid
	x = Spharmt(nlons,nlats,ntrunc,rsphere,gridtype='gaussian')
	lons,lats = np.meshgrid(x.lons, x.lats)


	if False:
		print("*"*80)
		print("SIN")
		print("*"*80)
		ug = np.sin(x.lats)	# Sinus generates significant errors
		ug.shape = (nlats,1)
		ug = ug*np.ones((nlats,nlons),dtype=np.float)
		vg = np.zeros((nlats,nlons),np.float)

		vrtspec, divspec =  x.getvrtdivspec(ug,vg)
		ugx, vgx = x.getuv(vrtspec, divspec)
		
		outputMinMaxSum(ug, "ug")
		outputMinMaxSum(ugx, "ugx")
		outputMinMaxSum(ug-ugx, "ug-ugx")

		c = 0
		for i in range(T+1):
			for j in range(i, T+1):
				print(vrtspec[c])
				c = c+1
			print("*"*80)
		sys.exit(1)



	print("*"*80)
	print("COS")
	print("*"*80)
	ug = np.cos(x.lats)	# Cosinus does not generate significant errors
	ug.shape = (nlats,1)
	ug = ug*np.ones((nlats,nlons),dtype=np.float)
	vg = np.zeros((nlats,nlons),np.float)

	vrtspec, divspec =  x.getvrtdivspec(ug,vg)
	ugx, vgx = x.getuv(vrtspec, divspec)
	
	outputMinMaxSum(ug, "ug")
	outputMinMaxSum(ugx, "ugx")
	outputMinMaxSum(ug-ugx, "ug-ugx")

