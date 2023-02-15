#! /usr/bin/env python3

import numpy as np
import sys
import math

from mule.postprocessing.SphereData import *


class SphereDataNormsPhysicalSpace:

    def __init__(self, filename_a = None, filename_b = None, i_output_file = None):

        if filename_b != None:
            self.compute_diff(filename_a, filename_b)

        if i_output_file is not None:
            self.write_file(i_output_pickle_file)



    def compute_diff(
            self,
            filename_a,
            filename_b,
        ):
        file_a = SphereData(filename_a, setup_physical=True)
        file_b = SphereData(filename_b, setup_physical=True)

        self.norm_l1_value = 0.0
        self.norm_l2_value = 0.0
        self.norm_linf_value = 0.0
        self.norm_rms_value = 0.0

        size_ref_j = file_a.data_physical.shape[0]
        size_ref_i = file_a.data_physical.shape[1]
        size_cmp_j = file_b.data_physical.shape[0]
        size_cmp_i = file_b.data_physical.shape[1]

        multiplier_j = (size_ref_j+1)/(size_cmp_j+1)
        multiplier_i = (size_ref_i+1)/(size_cmp_i+1)


        print ("Dimensions of reference solution: ", size_ref_i, size_ref_j)
        print ("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
        if not float(multiplier_i).is_integer() or not float(multiplier_j).is_integer() : 
            print ("Grids are not aligned")
            print ("Try to use (TODO) interpolation script")
            print ("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
            print ("Multipliers: ", multiplier_i, multiplier_j)
            raise Exception("Grids not properly aligned")

        multiplier_j = int(multiplier_j)
        multiplier_i = int(multiplier_i)

        print("Using multipliers (int): ", multiplier_i, multiplier_j)

        for j in range(0, size_cmp_j):
            for i in range(0, size_cmp_i):
                value = file_b.data_physical[j,i]-file_a.data_physical[j*multiplier_j,i*multiplier_i]

                # http://mathworld.wolfram.com/L1-Norm.html
                self.norm_l1_value += abs(value)
                # http://mathworld.wolfram.com/L2-Norm.html
                self.norm_l2_value += value*value
                # http://mathworld.wolfram.com/L-Infinity-Norm.html
                self.norm_linf_value = max(abs(value), self.norm_linf_value)

                # http://mathworld.wolfram.com/Root-Mean-Square.html
                self.norm_rms_value += value*value

        self.N = size_cmp_i*size_cmp_j

        # Compute sqrt() for Euklidian L2 norm
        self.norm_l2_value = math.sqrt(self.norm_l2_value)

        # RMS final sqrt(N) computation
        self.norm_rms_value  = math.sqrt(self.norm_rms_value/self.N)

        # resolution normalized L1 value
        self.res_norm_l1_value = self.norm_l1_value/float(self.N)


    def print(self, prefix=""):
        print(f"{prefix}norm l1: {self.norm_l1_value}")
        print(f"{prefix}norm l2: {self.norm_l2_value}")
        print(f"{prefix}norm linf: {self.norm_linf_value}")
        print(f"{prefix}norm rms: {self.norm_rms_value}")
        print(f"{prefix}res norm l1: {self.res_norm_l1_value}")



    def write_file(
            self,
            picklefile,
            tagprefix = None,
            verbosity = 0
        ):

        #
        # If picklefile is specified, write norm data to pickle file.
        # This can be later on further postprocessed!
        #
        import pickle

        if tagprefix != None:
            if tagprefix[-1] != ".":
                tagprefix += '.'
        else:
            tagprefix = ""

        pickle_data = {
            tagprefix+'N' : self.N,
            tagprefix+'norm_l1' : self.norm_l1_value,
            tagprefix+'norm_l2' : self.norm_l2_value,
            tagprefix+'norm_linf' : self.norm_linf_value,
            tagprefix+'norm_rms' : self.norm_rms_value,
        }

        # Write values for resolution neutral values into the .pickle files
        pickle_data.update({
            tagprefix+'res_norm_l1' : self.res_norm_l1_value,
            tagprefix+'res_norm_l2' : self.norm_rms_value,    # This is the RMS = L2
            tagprefix+'res_norm_linf' : self.norm_linf_value,    # No normalization required
            tagprefix+'res_norm_rms' : self.norm_rms_value,    # Already normalized
        })

        pickle_data['WARNING'] = "L1, L2 and RMS don't include scaling factors for different cell spacings around the sphere!!!"

        print(" + picklefile: "+str(picklefile))

        with open(picklefile, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(pickle_data, f)

        print("")
