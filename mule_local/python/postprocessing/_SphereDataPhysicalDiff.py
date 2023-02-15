#! /usr/bin/env python3

import numpy as np
import sys
import math

from mule.postprocessing.SphereData import *



class SphereDataPhysicalDiff:

    def __init__(self, filename_a = None, filename_b = None, i_output_file = None):

        if filename_b != None:
            self.compute_diff(filename_a, filename_b)

        if i_output_file is not None:
            self.write_file(i_output_pickle_file)



    def compute_diff(
            self,
            filename_a,
            filename_b,
            verbosity = 0
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


        if verbosity > 0:
            print("Dimensions of reference solution: ", size_ref_i, size_ref_j)
            print("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
        if not float(multiplier_i).is_integer() or not float(multiplier_j).is_integer() : 
            print("Grids are not aligned")
            print("Try to use (TODO) interpolation script")
            print("Dimensions of method under analysis: ", size_cmp_i, size_cmp_j)
            print("Multipliers: ", multiplier_i, multiplier_j)
            raise Exception("Grids not properly aligned")

        multiplier_j = int(multiplier_j)
        multiplier_i = int(multiplier_i)

        if verbosity > 0:
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
            tagname = None,
            verbosity = 0
        ):

        #
        # If picklefile is specified, write norm data to pickle file.
        # This can be later on further postprocessed!
        #
        import pickle

        if tagname != None:
            tagname += '.'
        else:
            tagname = ''

        pickle_data = {
            tagname+'N' : self.N,
            tagname+'diff.norm_l1' : self.norm_l1_value,
            tagname+'diff.norm_l2' : self.norm_l2_value,
            tagname+'diff.norm_linf' : self.norm_linf_value,
            tagname+'diff.norm_rms' : self.norm_rms_value,
            tagname+'diff.res_norm_l1' : self.res_norm_l1_value,
        }

        pickle_data['WARNING'] = "L1 and L2 don't include scaling factors for different cell spacings around the sphere!!!"

        print(" + picklefile: "+str(picklefile))

        with open(picklefile, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(pickle_data, f)

        print("")



if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("")
        print("Usage:")
        print("    "+sys.argv[0]+" [infile A] [infile B] [picklefile output file (optional)] [reference tagname]")
        print("")
        print("    infile A:")
        print("        First input .[sweet/csv] file with spectral or physical space data on the sphere")
        print("")
        print("    infile B:")
        print("        Second input .[sweet/csv] file with spectral or physical space data on the sphere")
        print("")
        print("    picklefile:")
        print("        If given, output is pickled into this file")
        print("        .diff.norm_l1")
        print("        .diff.norm_l2")
        print("        .diff.norm_linf")
        print("        .diff.norm_rms")
        print("        .diff.res_norm_l1")
        print("")
        print(" reference tagname:")
        print("        How to name value in .pickle file")
        print("")
        sys.exit(1)

    filename_a = sys.argv[1]
    filename_b = sys.argv[2]

    picklefile = None
    if len(sys.argv) > 3:
        picklefile = sys.argv[3]

    tagname = None
    if len(sys.argv) > 4:
        tagname = sys.argv[4]

    s = SphereDataPhysicalDiff()
    s.compute_diff(filename_a, filename_b)
    s.print()
    s.write_file(picklefile, tagname)
