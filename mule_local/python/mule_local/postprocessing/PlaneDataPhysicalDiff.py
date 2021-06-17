#! /usr/bin/env python3

import sys
import math

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy import signal

from mule_local.postprocessing.PlaneDataPhysical import *
from mule.InfoError import *



class PlaneDataPhysicalDiff(InfoError):

    def __init__(
        self,
        filename_reference = None,
        filename_b = None,
        params = ['interpolate']
    ):
        """
        Parameters:

        filname_reference: string
            Filename (full path) of reference data

        filename_b: string
            Filename (fullpath) of data to compare with


        params: list of strings

            'default':
                only support 1:1 resolution matching

            All the following methods assume cell-centered data.
            This is typically the case for the height data field.

            'reduce_ref_to_cmp':
                reduce reference solution to lower resolution by accumulating data from cells

            'interpolate':
                if resolutions don't match, interpolate to reference resolution (cubic interpolation)

            'interpolate_ref_to_cmp':
                Use spline interpolation to interpolate reference to comparison data

            'interpolate_cmp_to_ref':
                Use spline interpolation to interpolate comparison to reference data

            'spectral_ref_to_cmp':
                Use spectral (modal) projection
        """

        InfoError.__init__(self, "PlaneDataPhysicalDiff")
        self.params = params

        if filename_b != None:
            self.compute_diff(filename_reference, filename_b, params)

        pass



    def compute_diff(
            self,
            filename_a,
            filename_b,
            params,
        ):

        # Load first file_b, to avoid wasting time for filename_a if filename_b doesn't exist
        file_b = PlaneDataPhysical(filename_b)
        file_a = PlaneDataPhysical(filename_a)


        self.norm_l1_value = 0.0
        self.norm_l2_value = 0.0
        self.norm_linf_value = 0.0
        self.norm_rms_value = 0.0
        self.N = 0

        size_ref_j = len(file_a.data)
        size_ref_i = len(file_a.data[0])
        size_cmp_j = len(file_b.data)
        size_cmp_i = len(file_b.data[0])

        if size_ref_j == size_cmp_j and size_ref_i == size_cmp_i:# and False:
            self.info("Using 'default' method (matching resolutions)")

            data_ref = file_a.data
            data_cmp = file_b.data


        elif 'reduce_ref_to_cmp' in params:
            self.info("Using 'ref_reduce' method")

            #
            # Reduce reference solution, assuming that the
            # resolution is integer multiples of the one to compare with
            #

            data_ref = file_a.data
            data_cmp = file_b.data

            multiplier_j = size_ref_j/size_cmp_j
            multiplier_i = size_ref_i/size_cmp_i

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

            print("Using multipliers (int): ", multiplier_i, multiplier_j)

            data_ref_new = np.zeros(shape=data_cmp.shape)
            for j in range(0, size_cmp_j):
                for i in range(0, size_cmp_i):

                    data_ref_new[j,i] = 0.0
                    for sj in range(0, multiplier_j):
                        for si in range(0, multiplier_i):
                            data_ref_new[j,i] += data_ref[j*multiplier_j+sj, i*multiplier_i+si]

                    data_ref_new[j,i] /= multiplier_i*multiplier_j


            data_ref = data_ref_new


        elif ('interpolate' in params and size_ref_j > size_cmp_j and size_ref_i > size_cmp_i) or ('interpolate_cmp_to_ref' in params):
            self.info("Using 'interpolate_cmp_to_ref' method")

            #
            # Interpolate solution to resolution of reference solution.
            # Note, that this can also lead to a reduction in case of a lower resolution of the reference
            #

            # This is the code written by Pedro which uses spline interpolation

            # Comparison via interpolation
            # print("Interpolation")
            # A-grid REFERENCE (file1) - sweet outputs only A grids physical space

            data_ref = file_a.data
            data_cmp = file_b.data

            ny_ref = len(file_a.data)
            nx_ref = len(file_a.data[0])
            print("REF resolution: "+str(nx_ref)+", "+str(ny_ref))

            ny_cmp = len(file_b.data)
            nx_cmp = len(file_b.data[0])
            print("CMP resolution: "+str(nx_cmp)+", "+str(ny_cmp))

            dx_ref = 1.0/(nx_ref)
            dy_ref = 1.0/(ny_ref)
            x_ref = np.linspace(0, 1, nx_ref, endpoint=False)
            y_ref = np.linspace(0, 1, ny_ref, endpoint=False)
            x_ref += dx_ref/2    # Make it cell-centered
            y_ref += dy_ref/2
            X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

            # A-grid cmp file (file2)
            dx_cmp=1.0/nx_cmp
            dy_cmp=1.0/ny_cmp
            x_cmp = np.linspace(0, 1, nx_cmp, endpoint=False)
            y_cmp = np.linspace(0, 1, ny_cmp, endpoint=False)
            x_cmp += dx_cmp/2
            y_cmp += dy_cmp/2
            X_cmp, Y_cmp = np.meshgrid(x_cmp, y_cmp)

            # Create cubic interpolation of comparison
            # Data to be interpolated
            interp_spline = RectBivariateSpline(x_cmp, y_cmp, data_cmp)

            # Compute reduced reference resolution
            # Target grid
            data_cmp_high = interp_spline(x_ref, y_ref)

            data_cmp = data_cmp_high
            data_ref = data_ref


        elif ('interpolate' in params and size_ref_j < size_cmp_j and size_ref_i < size_cmp_i) or ('interpolate_ref_to_cmp' in params):
            self.info("Using 'interpolate_ref_to_cmp' method")

            #
            # Interpolate: REF => resolution(CMP)
            #

            #
            # Interpolate solution to resolution of reference solution.
            # Note, that this can also lead to a reduction in case of a lower resolution of the reference
            #

            # This is the code written by Pedro which uses spline interpolation

            # Comparison via interpolation
            # print("Interpolation")
            # A-grid REFERENCE (file1) - sweet outputs only A grids physical space

            data_ref = file_a.data
            data_cmp = file_b.data

            ny_ref = len(file_a.data)
            nx_ref = len(file_a.data[0])
            print("REF resolution: "+str(nx_ref)+", "+str(ny_ref))

            ny_cmp = len(file_b.data)
            nx_cmp = len(file_b.data[0])
            print("CMP resolution: "+str(nx_cmp)+", "+str(ny_cmp))

            dx_ref = 1.0/nx_ref
            dy_ref = 1.0/ny_ref
            x_ref = np.linspace(0, 1, nx_ref, endpoint=False)
            y_ref = np.linspace(0, 1, ny_ref, endpoint=False)
            x_ref += dx_ref/2    # Make it cell-centered, this significantly reduces the errors
            y_ref += dy_ref/2
            X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

            # A-grid cmp file (file2)
            dx_cmp=1.0/nx_cmp
            dy_cmp=1.0/ny_cmp
            x_cmp = np.linspace(0, 1, nx_cmp, endpoint=False)
            y_cmp = np.linspace(0, 1, ny_cmp, endpoint=False)
            x_cmp += dx_cmp/2
            y_cmp += dy_cmp/2
            X_cmp, Y_cmp = np.meshgrid(x_cmp, y_cmp)

            # Create cubic interpolation of reference file
            interp_spline = RectBivariateSpline(x_ref, y_ref, data_ref)

            # Compute reduced reference resolution
            data_ref_low = interp_spline(x_cmp, y_cmp)

            data_ref = data_ref_low


        #elif 'spectral_ref_to_cmp' in params:
        elif 'spectral_ref_to_cmp' in params or 'spectral' in params:

            self.info("Using 'spectral_ref_to_cmp' method")

            # This is the code written by Pedro which uses spline interpolation

            # Comparison via interpolation
            # print("Interpolation")
            # A-grid REFERENCE (file1) - sweet outputs only A grids physical space

            data_ref = file_a.data
            data_cmp = file_b.data

            """
            #if True:
            if False:
                data_ref = np.zeros((8, 8))
                data_cmp = np.zeros((6, 6))
                for i in range(len(data_ref[0])):
                    for j in range(len(data_ref)):
                        data_ref[j,i] = 1.0

                        a = 1.0
                        a *= math.sin(2.0*i*math.pi*2.0/len(data_ref[0]))
                        a *= math.cos(2.0*j*math.pi*2.0/len(data_ref))
                        data_ref[j,i] += a

                        a = 1.0
                        a *= math.sin(i*math.pi*2.0/len(data_ref[0]))
                        a *= math.cos(j*math.pi*2.0/len(data_ref))
                        data_ref[j,i] += a

                for i in range(len(data_cmp[0])):
                    for j in range(len(data_cmp)):
                        data_cmp[j,i] = 1.0

                        a = 1.0
                        a *= math.sin(2.0*i*math.pi*2.0/len(data_cmp[0]))
                        a *= math.cos(2.0*j*math.pi*2.0/len(data_cmp))
                        data_cmp[j,i] += a

                        a = 1.0
                        a *= math.sin(i*math.pi*2.0/len(data_cmp[0]))
                        a *= math.cos(j*math.pi*2.0/len(data_cmp))
                        data_cmp[j,i] += a
            """

            ny_ref = len(data_ref)
            nx_ref = len(data_ref[0])
            res_ref = nx_ref*ny_ref
            print("REF resolution: "+str(nx_ref)+", "+str(ny_ref))

            ny_cmp = len(data_cmp)
            nx_cmp = len(data_cmp[0])
            res_cmp = nx_cmp*ny_cmp
            print("CMP resolution: "+str(nx_cmp)+", "+str(ny_cmp))

            if nx_cmp & 1 or ny_cmp & 1:
                raise Exception("Only even resolutions allowed")

            #
            # Spectral rescaling
            #

            shift_ref_i = -1.0/size_ref_i*0.5
            shift_ref_j = -1.0/size_ref_j*0.5

            shift_cmp_j = 1.0/size_cmp_j*0.5
            shift_cmp_i = 1.0/size_cmp_i*0.5

            #if False:
            if True:
                #
                # RFFT
                #

                # Convert reference data to spectrum
                data_ref_spec = np.fft.rfft2(data_ref)

                # Dummy transformation to get matching spectral size
                data_ref_new_spec = np.fft.rfft2(np.zeros(shape=data_cmp.shape))


                specx = data_ref_spec.shape[1]
                specy = data_ref_spec.shape[0]

                newspecx = data_ref_new_spec.shape[1]
                newspecy = data_ref_new_spec.shape[0]

                data_ref_new_spec[0:newspecy//2, 0:newspecx] = data_ref_spec[0:newspecy//2, 0:newspecx]

                d = specy - newspecy//2
                nd = newspecy - newspecy//2
                data_ref_new_spec[nd:nd+newspecy//2, 0:newspecx] = data_ref_spec[d:d+newspecy//2, 0:newspecx]

                """
                print("*"*80)
                print(data_ref_spec)
                print("*"*80)
                print(data_ref_new_spec)
                print("*"*80)
                """

                def shift(
                        spec_data,    # spectral data
                        sh_x,        # shift in x direction within domain [0;1]
                        sh_y,        # shift in y direction within domain [0;1]
                    ):

                    # return data
                    ret_data = np.zeros(spec_data.shape, dtype=complex)

                    for iy in range(spec_data.shape[0]):
                        # compute mode
                        if iy < spec_data.shape[0]//2:
                            ky = iy
                        else:
                            ky = iy-spec_data.shape[0]

                        for ix in range(spec_data.shape[1]):
                            # compute mode
                            kx = ix

                            ret_data[iy,ix] = spec_data[iy,ix]
                            #ret_data[iy,ix] *= np.exp(1j*2.0*math.pi*kx*sh_x)
                            #ret_data[iy,ix] *= np.exp(1j*2.0*math.pi*ky*sh_y)
                            ret_data[iy,ix] *= np.exp(1j*2.0*math.pi*(kx*sh_x + ky*sh_y))

                    return ret_data

                # Shift high res reference data to 0,0
                data_ref_new_spec = shift(data_ref_new_spec, shift_ref_i+shift_cmp_i, shift_ref_j+shift_cmp_j)

                data_ref_new = np.fft.irfft2(data_ref_new_spec)
                data_ref_new *= res_cmp/res_ref


            else:
                #
                # FFT
                #

                # Convert reference data to spectrum
                data_ref_spec = np.fft.fft2(data_ref)

                # Dummy transformation to get matching spectral size
                data_ref_new_spec = np.fft.fft2(np.zeros(shape=data_cmp.shape))



                specx = data_ref_spec.shape[1]
                specy = data_ref_spec.shape[0]

                newspecx = data_ref_new_spec.shape[1]
                newspecy = data_ref_new_spec.shape[0]

                dy = specy - newspecy//2
                ndy = newspecy - newspecy//2

                dx = specx - newspecx//2
                ndx = newspecx - newspecx//2

                data_ref_new_spec[    0:newspecy//2,            0:newspecx//2    ] = data_ref_spec[    0:newspecy//2,        0:newspecx//2    ]
                data_ref_new_spec[    ndy:ndy+newspecy//2,    0:newspecx//2    ] = data_ref_spec[    dy:dy+newspecy//2,    0:newspecx//2    ]

                data_ref_new_spec[    0:newspecy//2,            ndx:ndx+newspecx//2    ] = data_ref_spec[    0:newspecy//2,        dx:dx+newspecx//2    ]
                data_ref_new_spec[    ndy:ndy+newspecy//2,    ndx:ndx+newspecx//2    ] = data_ref_spec[    dy:dy+newspecy//2,    dx:dx+newspecx//2    ]

                if False:
                    print("*"*80)
                    print(data_ref_spec)
                    print("*"*80)
                    print(data_ref_new_spec)
                    print("*"*80)


                def shift(
                        spec_data,    # spectral data
                        sh_x,        # shift in x direction within domain [0;1]
                        sh_y,        # shift in y direction within domain [0;1]
                    ):

                    # return data
                    ret_data = np.zeros(spec_data.shape, dtype=complex)

                    for iy in range(spec_data.shape[0]):
                        # compute mode
                        if iy < spec_data.shape[0]//2:
                            ky = iy
                        else:
                            ky = iy-spec_data.shape[0]

                        for ix in range(spec_data.shape[1]):
                            # compute mode
                            if ix < spec_data.shape[1]//2:
                                kx = ix
                            else:
                                kx = ix-spec_data.shape[1]

                            ret_data[iy,ix] = spec_data[iy,ix]
                            ret_data[iy,ix] *= np.exp(1j*2.0*math.pi*kx*sh_x)
                            ret_data[iy,ix] *= np.exp(1j*2.0*math.pi*ky*sh_y)

                    return ret_data

                # Shift high res reference data to 0,0
                data_ref_new_spec = shift(data_ref_new_spec, shift_ref_i+shift_cmp_i, shift_ref_j+shift_cmp_j)

                data_ref_new = np.fft.ifft2(data_ref_new_spec)
                data_ref_new *= res_cmp/res_ref

                # Restrict to real data
                data_ref_new = np.real(data_ref_new)


            # The axes might be wrong.
            #
            # Swap axis=1 and axis=0 if there are some errors.
            #
            #f = signal.resample(data_ref, nx_cmp, t=None, axis=0)
            #data_ref_new = signal.resample(f, ny_cmp, t=None, axis=1)

            data_ref = data_ref_new


        else:
            print("")
            print("No supported method provided in '"+(','.join(params))+"'")
            print("")
            raise Exception("No supported method provided in '"+(','.join(params))+"'")


        #
        # Grids have same resolution
        #
        for j in range(0, data_cmp.shape[0]):
            for i in range(0, data_cmp.shape[1]):
                value = data_cmp[j,i] - data_ref[j,i]

                # http://mathworld.wolfram.com/L1-Norm.html
                self.norm_l1_value += abs(value)

                # http://mathworld.wolfram.com/L2-Norm.html
                self.norm_l2_value += value*value

                # http://mathworld.wolfram.com/L-Infinity-Norm.html
                self.norm_linf_value = max(abs(value), self.norm_linf_value)

                # http://mathworld.wolfram.com/Root-Mean-Square.html
                self.norm_rms_value += value*value

                self.N += 1



        # Compute sqrt() for Euklidian L2 norm
        self.norm_l2_value = math.sqrt(self.norm_l2_value)

        # RMS final sqrt(N) computation
        self.norm_rms_value  = math.sqrt(self.norm_rms_value/float(self.N))

        # resolution normalized L1 value
        self.res_norm_l1_value = self.norm_l1_value/float(self.N)




    def print(self):
        print("")
        print(" + norm l1: "+str(self.norm_l1_value))
        print(" + norm l2: "+str(self.norm_l2_value))
        print(" + norm linf: "+str(self.norm_linf_value))
        print(" + norm rms: "+str(self.norm_rms_value))
        print(" + res norm l1: "+str(self.res_norm_l1_value))
        print("")



    def write_file(
            self,
            picklefile,
            tagname = None
        ):

        #
        # If picklefile is specified, write norm data to pickle file.
        # This can be later on further postprocessed!
        #
        if picklefile != None:
            import pickle

            if tagname != None:
                tagname += '.'
            else:
                tagname = ''

            pickle_data = {
                tagname+'N' : self.N,
                tagname+'norm_l1' : self.norm_l1_value,
                tagname+'norm_l2' : self.norm_l2_value,
                tagname+'norm_linf' : self.norm_linf_value,
                tagname+'norm_rms' : self.norm_rms_value,
            }

            # Write values for resolution neutral values into the .pickle files
            pickle_data.update({
                tagname+'res_norm_l1' : self.res_norm_l1_value,
                tagname+'res_norm_l2' : self.norm_rms_value,    # This is the RMS = L2
                tagname+'res_norm_linf' : self.norm_linf_value,    # No normalization required
                tagname+'res_norm_rms' : self.norm_rms_value,    # Already normalized
            })

            print(" + picklefile: "+str(picklefile))

            with open(picklefile, 'wb') as f:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle.dump(pickle_data, f)

        print("")



if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("")
        print("Usage:")
        print("    "+sys.argv[0]+" [infile A] [infile B] [parameters] [picklefile output file (optional)] [reference tagname (optional)]")
        print("")
        print("    infile A:")
        print("        First input .csv file with physical space data on the plane")
        print("")
        print("    infile B:")
        print("        Second input .csv file with physical space data on the plane")
        print("")
        print("    parameters:")
        print("        'interpolate':")
        print("            use interpolation if resolutions don't match")
        print("")
        print("    picklefile:")
        print("        If given, output is pickled into this file")
        print("        diff.error_l1")
        print("        diff.error_l2")
        print("        diff.error_linf")
        print("        diff.error_rms")
        print("")
        print(" reference tagname:")
        print("        How to name value in .pickle file")
        print("")
        sys.exit(1)



    filename_a = sys.argv[1]
    filename_b = sys.argv[2]

    params = []
    if len(sys.argv) > 3:
        params = sys.argv[3].split(' ')

    picklefile = None
    if len(sys.argv) > 4:
        picklefile = sys.argv[4]

    tagname = None
    if len(sys.argv) > 5:
        tagname = sys.argv[5]

    print("")
    print("parameters: "+str(params))
    print("picklefile: "+str(picklefile))
    print("tagname: "+str(tagname))
    print("")

    s = PlaneDataPhysicalDiff()
    s.compute_diff(filename_a, filename_b, params)
    s.print()

    if picklefile != None:
        s.write_file(picklefile, tagname)
