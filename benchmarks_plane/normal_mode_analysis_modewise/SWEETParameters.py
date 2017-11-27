#! /usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import os
import sys
import stat
import math
from SWEETJobGeneration import *



class SWEETCompilePlane:

	def __init__(self):
		p = SWEETJobGeneration()
		print(p.compile)
		p.compile = 'plane'

		p.compile.program = 'swe_plane'
		
		p.compile.plane_or_sphere = 'plane'
		p.compile.plane_spectral_space = 'enable'
		p.compile.plane_spectral_dealiasing = 'enable'
		p.compile.sphere_spectral_space = 'disable'
		p.compile.sphere_spectral_dealiasing = 'disable'


		p.compile.compiler = 'gnu'
		self = p



class SWEETParamsEarth:

	def __init__(self):
		p = SWEETJobGeneration()

		p.compile = 'plane'

		p.compile.program = 'swe_plane'
		
		p.compile.plane_or_sphere = 'plane'
		p.compile.plane_spectral_space = 'enable'
		p.compile.plane_spectral_dealiasing = 'enable'
		p.compile.sphere_spectral_space = 'disable'
		p.compile.sphere_spectral_dealiasing = 'disable'


		p.compile.compiler = 'gnu'
		self = p

