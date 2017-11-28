#! /usr/bin/env python3
#
#
#   Pre-setup specific test parameters
#
#
#
#
#
#---------------------------------------------------

import matplotlib
matplotlib.use('agg')

import os
import sys
import stat
import math

def CompileSWEPlane(p):
	p.compile.program = 'swe_plane'
	p.compile.plane_or_sphere = 'plane'
	p.compile.plane_spectral_space = 'enable'
	p.compile.plane_spectral_dealiasing = 'enable'
	p.compile.sphere_spectral_space = 'disable'
	p.compile.sphere_spectral_dealiasing = 'disable'
	p.compile.compiler = 'gnu'
	return p	



def RuntimeSWEPlaneNondimParameters(p):
	p.runtime.g = 1
	p.runtime.f = 1
	p.runtime.h = 1
	p.runtime.domain_size = 1
	return p	

def RuntimeSWEPlaneEarthParam(p):
	p.runtime.g = 9.80616
	p.runtime.f = 0.00014584
	p.runtime.h = 10000
	p.runtime.domain_size = 40031555.8928087
	return p	

def EnableGUI(p):
	p.runtime.gui = 'enable'
	p.compile.gui = 'enable'
	return p	


def DisableGUI(p):
	p.runtime.gui = 'disable'
	p.compile.gui = 'disable'
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
