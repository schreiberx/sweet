#! /usr/bin/env python3

import numpy as np
import sys
import math

from mule.postprocessing.ScalarData import *


class ScalarDataNormsPhysicalSpace:

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
        file_a = ScalarData(filename_a)
        file_b = ScalarData(filename_b)

        self.norm_l1_value = 0.0
        self.norm_l2_value = 0.0
        self.norm_linf_value = 0.0
        self.norm_rms_value = 0.0

        size_ref = file_a.data.size
        size_cmp = file_b.data.size

        if not (size_ref == size_cmp):
            raise Exception("Files do not have the same data size")


        print ("Dimension of solution: ", size_cmp)
        for i in range(0, size_cmp):
            value = file_b.data[i]-file_a.data[i]

            # http://mathworld.wolfram.com/L1-Norm.html
            self.norm_l1_value += abs(value)
            # http://mathworld.wolfram.com/L2-Norm.html
            self.norm_l2_value += value*value
            # http://mathworld.wolfram.com/L-Infinity-Norm.html
            self.norm_linf_value = max(abs(value), self.norm_linf_value)

            # http://mathworld.wolfram.com/Root-Mean-Square.html
            self.norm_rms_value += value*value

        self.N = size_cmp

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

        print(" + picklefile: "+str(picklefile))

        with open(picklefile, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(pickle_data, f)

        print("")
