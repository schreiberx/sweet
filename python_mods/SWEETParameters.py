#! /usr/bin/env python3
#
#
#   Pre-setup specific test parameters
#
#
#   Pedro Peixoto <pedrosp@ime.usp.br>
#
#
#---------------------------------------------------

import os
import sys
import stat
import math


class EarthMKSDimensions:
	name = "EarthMKSDimensions"
	day = 86400
	g = 9.80616
	f = 0.00014584	#2 omega
	omega = 0.00007292
	erad = 6371220


def CompileSWEPlane(p):
	p.compile.program = 'swe_plane'
	p.compile.plane_or_sphere = 'plane'
	p.compile.plane_spectral_space = 'enable'
	#p.compile.plane_spectral_dealiasing = 'enable'
	p.compile.plane_spectral_dealiasing = 'disable'
	p.compile.sphere_spectral_space = 'disable'
	p.compile.sphere_spectral_dealiasing = 'disable'
	p.compile.compiler = 'gnu'
	p.compile.threading = 'omp'
	p.compile.quadmath = 'disable'
	return p


def RuntimeSWEPlaneNondimParameters(p):
	p.runtime.g = 1
	p.runtime.f = 1
	p.runtime.h = 1
	p.runtime.domain_size = 1
	return p

def RuntimeSWEPlaneEarthParam(p):
	s = EarthMKSDimensions()
	p.runtime.g = s.g
	p.runtime.f = s.f
	p.runtime.h = 10000
	p.runtime.domain_size = 2.0*math.pi*s.erad # 40031555.8928087
	return p


def EnableGUI(p):
	p.runtime.gui = 'enable'
	p.compile.gui = 'enable'
	return p


def DisableGUI(p):
	p.runtime.gui = 'disable'
	p.compile.gui = 'disable'
	return p


def SetupFDCMethods(p):
	p.runtime.staggering = 1
	p.runtime.spectralderiv = 0

	#p.compile.plane_spectral_space = 'disable'
	p.compile.plane_spectral_dealiasing = 'disable'
	p.compile.libfft = 'enable'
	return p


def SetupSpectralMethods(p):
	p.runtime.staggering = 0
	p.runtime.spectralderiv = 1

	p.compile.plane_spectral_space = 'enable'
	p.compile.plane_spectral_dealiasing = 'enable'
	return p



class Unstablejet:
	name = "Unstablejet"
	s = EarthMKSDimensions()
	#Domain
	x_min = 0
	x_max = 2.0*math.pi*s.erad/1000/1000 #1000km
	y_min = 0
	y_max = 2.0*math.pi*s.erad/1000/1000 #1000km

class ParameterFilename:
	def __init__(self, filename):
		self.filename = filename
		self.basename = os.path.basename(filename)
		self.name, self.ext = os.path.splitext(self.basename)
		self.dirname = os.path.dirname(filename)
		self.extract_all()
		self.heading = "Variable \t Method \t MethodPaper \t Time \t dt \t"
		self.details = self.outname+"\t"+self.method+"\t"+self.method_paper+"\t"+self.time+"\t"+self.timestep+"\t"
		self.prefix = self.outname+"_"+self.method+"_t"+self.time+"_dt"+self.timestep
	def extract_all(self):
		self.setup()
		self.extract_var()
		self.extract_method()
		self.extract_time()

	def setup(self):
		#Check if plotting in standard naming convention
		if 'script_' in self.filename:
			self.stdpat = 1
		else:
			self.stdpat = 0

	def extract_var(self):
		#Check variable to be plotted
		if 'diag_vort' in self.basename:
			self.title = "Vorticity Deviation \n"
			self.outname = "Vort "
			self.var = 'vort'

		elif 'prog_h' in self.basename:
			self.title = "Depth Deviation (m) \n"
			self.outname = "Depth "
			self.var = "depth"

		else:
			pos1 = self.basename.find('output_')
			pos2 = self.basename.find('_t')
			if pos1 < 0:
				self.var = self.basename
			else:
				self.var = self.filename[pos1+8:pos2]
			self.title = self.var
			self.outname = serf.var


	def extract_method(self):
		#Method
		if self.stdpat == 1:
			pos1 = self.filename.find('_tsm_')
			pos2 = self.filename.find('_tso')
			self.method_original = self.filename[pos1+5:pos2]
			self.method_paper = self.method_original
			if self.method_original == "l_cn_na_sl_nd_settls":
				self.method_paper = "SL-SI-SETTLS"
			elif self.method_original == "l_rexi_na_sl_nd_settls":
				self.method_paper = "SL-EXP-SETTLS"
			elif self.method_original == "l_rexi_na_sl_nd_etdrk":
				self.method_paper = "SL-ETD2RK"
			elif self.method_original == "l_rexi_n_etdrk":
				self.method_paper = "ETD2RK"
			elif self.method_original == "ln_erk":
				if 'ref' in self.filename:
					self.method_paper = "REF"
				else:
					self.method_paper = "RK-FDC"
		else:
			self.method_original = self.dirname 
			self.method_paper = self.dirname
		self.method = self.method_original

	def extract_time(self):
		if self.stdpat == 1:
			#Time in days
			pos1 = self.filename.find('output')
			name = self.filename[pos1:]
			pos2 = name.find('_t')
			pos3 = self.filename.find('.csv')
			time1 = self.filename[pos1+pos2+2:pos3]
			time1 = float(time1)
			self.time = str(round(time1 / 86400, 2))

			
			#Time step in sec
			pos1 = self.filename.find('_C')
			pos2 = self.filename.find('_R')
			self.timestep=self.filename[pos1+2:pos2]

		else:
			self.time = ""
			self.timestep = ""

	

		
