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
	p.compile.plane_spectral_dealiasing = 'enable'
	p.compile.sphere_spectral_space = 'disable'
	p.compile.sphere_spectral_dealiasing = 'disable'
	p.compile.compiler = 'gnu'
	p.compile.threading = 'omp'
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








# ---old code ----

#from SWEETJobGeneration import *

#class SWEETSpecificTest(SWEETJobGeneration):
#
#	def __init__(self): 
#		self.p = SWEETJobGeneration() #For some odd reason, I need self.p, not just self :-(
#		self.p.compile.program="test"



#class SWEETEditTests(SWEETSpecificTest):

#	def __init__(self, p): 
#		self.p.compile.program="test2"
#		self = p
