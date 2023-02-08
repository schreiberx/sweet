#! /usr/bin/env python3

#
# Author: Martin Schreiber
# Email: M.Schreiber@exeter.ac.uk
# Date: 2017-06-16
#

import sys
import os
import copy

d = os.path.dirname(os.path.realpath(__file__))
sys.path.append(d+'/..')
import EFloat as ef
sys.path.pop()

import inspect



class TREXI_GaussianCoefficients:

    def __init__(self, floatmode = None):
        self.reset(floatmode = floatmode)

    def setupFloatmode(self):
        self.efloat = ef.EFloat(self.floatmode, self.mpfloat_digits)
        self.floatmode = self.efloat.floatmode


    def reset(self, floatmode = None):
        self.mpfloat_digits = None

        self.floatmode = floatmode

        # Reset floatmode in case that this was requested via a program parameter
        self.setupFloatmode()

        # Number of support points / Number of basis functions
        self.N = None

        # Basis function distance
        self.basis_function_spacing = None

        # Shift for rational approximated basis function
        self.basis_function_rat_shift = None


        self.weights = []
        self.weights_str = []
        self.weights_cplx = []

        # support points
        self.x0s = []



    def load_orig_ratgaussian_poles(self):
        self.reset()

        self.basis_function_rat_shift = self.efloat.to('-4.31532151087502402475593044073320925235748291015625')
        self.basis_function_spacing = self.efloat.to(1.0)

        self.N = 23

        weights_tmp = [
            ['-1.0845749544592896e-7',    '2.77075431662228e-8'],
            ['1.858753344202957e-8',    '-9.105375434750162e-7'],
            ['3.6743713227243024e-6',    '7.073284346322969e-7'],
            ['-2.7990058083347696e-6',    '0.0000112564827639346'],
            ['0.000014918577548849352',    '-0.0000316278486761932'],
            ['-0.0010751767283285608',    '-0.00047282220513073084'],
            ['0.003816465653840016',    '0.017839810396560574'],
            ['0.12124105653274578',        '-0.12327042473830248'],
            ['-0.9774980792734348',        '-0.1877130220537587'],
            ['1.3432866123333178',        '3.2034715228495942'],
            ['4.072408546157305',        '-6.123755543580666'],
            ['-9.442699917778205',        '0.'],
            ['4.072408620272648',        '6.123755841848161'],
            ['1.3432860877712938',        '-3.2034712658530275'],
            ['-0.9774985292598916',        '0.18771238018072134'],
            ['0.1212417070363373',        '0.12326987628935386'],
            ['0.0038169724770333343',    '-0.017839242222443888'],
            ['-0.0010756025812659208',    '0.0004731874917343858'],
            ['0.000014713754789095218',    '0.000031358475831136815 '],
            ['-2.659323898804944e-6',    '-0.000011341571201752273'],
            ['3.6970377676364553e-6',    '-6.517457477594937e-7'],
            ['3.883933649142257e-9',    '9.128496023863376e-7'],
            ['-1.0816457995911385e-7',    '-2.954309729192276e-8'],
        ]
        self.weights_cplx= []
        for w in weights_tmp:
            self.weights_cplx.append(self.efloat.cplx(self.efloat.to(w[0]), self.efloat.to(w[1])))

        # NOTE! We use the complex conjugates here which reflects the formulation
        # of the support points of the rational functions in this implementation
        self.weights_cplx = list(map(self.efloat.conj, self.weights_cplx))

        self.x0s = self.compute_basis_support_points()




    #
    # return the support points of the basis functions
    #
    def compute_basis_support_points(self):
        if self.N & 1 == 0:
            return [(self.efloat.to(-int(self.N/2)) + (self.efloat.to(i)+self.efloat.to(0.5)))*self.basis_function_spacing for i in range(self.N)]
        else:
            return [(self.efloat.to(-int(self.N/2)) + self.efloat.to(i))*self.basis_function_spacing for i in range(self.N)]
