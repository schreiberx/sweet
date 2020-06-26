#! /usr/bin/env python3

import sys, os
import numpy as np
import matplotlib.pyplot as plt

from mule.postprocessing.JobData import JobData
from mule_local.postprocessing.SphereDataSpectral import SphereDataSpectral
import mule_local.postprocessing.shtnsfiledata as shtnsfiledata


debug_active = False

class postprocessing_swe:
    def __init__(self):
        self.job_data = None
        self.job_data_flattened = None
        self.rsphere = None
        self.grav = None

        self.phi_pert_spec = None
        self.vrt_spec = None
        self.div_spec = None

        self.sh = None
        self.sh_aa = None


    def setup(self,
             i_job_dir,
             i_phi_file,
             i_vrt_file,
             i_div_file
    ):
        self.job_dir = i_job_dir
        self.phi_pert_file = i_phi_file
        self.vrt_file = i_vrt_file
        self.div_file = i_div_file

        # Load job's data
        self.job_data = JobData(self.job_dir)
        self.job_data_flattened = self.job_data.get_flattened_data()

        # TODO: Determine this automagically from job data
        self.rsphere = 6371220
        self.grav = 9.80616
        self.coriolis = 7.292*1e-5
        self.h0 = 10000

        # Load fields
        self.phi_pert_spec = SphereDataSpectral(self.phi_pert_file, setup_physical=False)
        self.vrt_spec = SphereDataSpectral(self.vrt_file, setup_physical=False)
        self.div_spec = SphereDataSpectral(self.div_file, setup_physical=False)

        # Setup transformations without anti-aliasing (lower physical resolution)
        self.sh = shtnsfiledata.shtnsfiledata(rsphere=self.rsphere)
        self.sh.setup(self.phi_pert_spec.file_info, anti_aliasing=False)

        # Setup transformations *with* anti-aliasing
        self.sh_aa = shtnsfiledata.shtnsfiledata(rsphere=self.rsphere)
        self.sh_aa.setup(self.phi_pert_spec.file_info, anti_aliasing=True)

        u_phys_data, v_phys_data = self.sh_aa.vrtdiv2uv(self.vrt_spec.data_spectral, self.div_spec.data_spectral)

    def plot_physical_field_data_only(
        self,
        i_tag,
        i_data_phys,
        i_output_filename
    ):
        numrows = i_data_phys.shape[0]
        numcols = i_data_phys.shape[1]

        #
        # Use a particular extension for the imshow
        # The rows need to be reversed for the contour plot for whatever reason
        #
        plt_extent = (-0.5, numcols-0.5, -0.5, numrows-0.5)

        plt.imshow(i_data_phys, extent=plt_extent)

        if i_tag == "h_pert":
            e = 50

        elif i_tag == "phi_pert":
            e = 50

        elif i_tag == "vrt":
            """
            Contours from
            Galewsky, J., Scott, R. K., & Polvani, L. M. (2004). An initial-value problem for testing numerical models of the global shallow-water equations. Tellus, Series A: Dynamic Meteorology and Oceanography, 56(5), 429–440. https://doi.org/10.1111/j.1600-0870.2004.00071.x
            Figure 1, Page 3
            
            They use the potential vorticity, h*vrt and contours e=0.2*\Omega/H
            Using just the vorticity, we get e=0.2*\Omega 
            """
            e = 0.2*self.coriolis

        elif i_tag == "div":
            e = 1e-7

        else:
            raise Exception("Unknown tag "+i_tag)

        num_contours = 30

        #
        # Interpolation zoom factor
        #
        contour_interpolation_zoom_factor = 1

        if contour_interpolation_zoom_factor != 1:
            from scipy import ndimage
            contour_data = ndimage.zoom(i_data_phys, mode='nearest', order=3, zoom=contour_interpolation_zoom_factor)
        else:
            contour_data = i_data_phys

        # For some reason, the rows need to be flipped for this plot
        plt_extent = (-0.5, numcols-0.5, numrows-0.5, -0.5)

        levels = np.linspace(e, e * num_contours, num_contours, endpoint=True)
        print("positive contour levels: "+str(levels))
        plt.contour(contour_data, levels=levels, extent=plt_extent, linestyles='solid', linewidths=0.2, colors='black')

        levels = np.linspace(-e * num_contours, -e, num_contours, endpoint=True)
        print("negative contour levels: "+str(levels))
        plt.contour(contour_data, levels=levels, extent=plt_extent, linestyles='dashed', linewidths=0.2, colors='black')

    def plot_physical_field(
        self,
        i_tag,
        i_data,
        i_output_filename,
        i_title
    ):
        print("Plotting ", i_output_filename)
        plt.close()

        self.plot_physical_field_data_only(i_tag, i_data, i_output_filename)
        plt.title(i_title)

        plt.tight_layout()
        plt.savefig(i_output_filename)

    def plot_physical_fields(self):
        # Compute u and v without anti-aliasing
        u_phys_data, v_phys_data = self.sh.vrtdiv2uv(self.vrt_spec.data_spectral, self.div_spec.data_spectral)

        vrt_phys_data = self.sh.spec2phys(self.vrt_spec.data_spectral)
        div_phys_data = self.sh.spec2phys(self.div_spec.data_spectral)

        # Compute u and v without anti-aliasing
        gh_pert_phys_data = self.sh.spec2phys(self.phi_pert_spec.data_spectral)

        # Compute height
        h_pert_phys_data = gh_pert_phys_data/self.grav


        ##############################
        # Plot the height field
        ##############################

        # Doesn't really exist, but just call it somehow
        input_filename = self.phi_pert_file.replace("phi", "h")

        output_filename = input_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')

        title = input_filename.replace('.pdf', '')
        _, title = os.path.split(title)

        self.plot_physical_field(
            "h_pert",
            h_pert_phys_data,
            output_filename,
            title
        )


        ##############################
        # Plot the vorticity field
        ##############################

        input_filename = self.vrt_file

        output_filename = input_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')

        title = input_filename.replace('.pdf', '')
        _, title = os.path.split(title)

        self.plot_physical_field(
            "vrt",
            vrt_phys_data,
            output_filename,
            title
        )


        ##############################
        # Plot the divergence field
        ##############################

        input_filename = self.div_file

        output_filename = input_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')

        title = input_filename.replace('.pdf', '')
        _, title = os.path.split(title)

        self.plot_physical_field(
            "div",
            div_phys_data,
            output_filename,
            title
        )

    def _ke_spectrum_dist_bucket(self, real_m, verbose=False):
        """
        Compute the bucket and distribution of mode real_m which is not integer

        Now we need to split things up into buckets, 0th mode bucket

        m = 0 ... 0.5 ... 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ...
            |      |               |               |
            | bck0 |    bucket1    |    bucket2    |   bucket3
            |      |               |               |

        Change the real_m to:

        m = 0.5 .. 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ... 3.5 ...
             |      |               |               |
             | bck0 |    bucket1    |    bucket2    |   bucket3
             |      |               |               |
        """

        if verbose:
            print(" +++ using real mode ", real_m)

        real_mh = real_m + 0.5

        # Now things are easy...

        bucket_a_num = int(real_mh - 0.5)
        bucket_b_num = int(real_mh + 0.5)

        bucket_a_weight = 1.0 - (real_mh - 0.5 - bucket_a_num)
        bucket_b_weight = real_mh + 0.5 - bucket_b_num

        if verbose:
            print(" +++ bucket ", bucket_a_num, " gets ", bucket_a_weight)
            print(" +++ bucket ", bucket_b_num, " gets ", bucket_b_weight)

        assert(bucket_a_num >= 0)
        assert(bucket_b_num == bucket_a_num + 1)
        assert(bucket_a_weight >= 0 and bucket_a_weight <= 1.0)
        assert np.allclose(bucket_a_weight + bucket_b_weight, 1.0)
        return bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight

    def _ke_spectrum_dist_bucket_array(self, real_m, verbose=False):
        """
        Compute the bucket and distribution of mode real_m which is not integer

        Now we need to split things up into buckets, 0th mode bucket

        m = 0 ... 0.5 ... 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ...
            |      |               |               |
            | bck0 |    bucket1    |    bucket2    |   bucket3
            |      |               |               |

        Change the real_m to:

        m = 0.5 .. 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ... 3.5 ...
             |      |               |               |
             | bck0 |    bucket1    |    bucket2    |   bucket3
             |      |               |               |
        """

        if verbose:
            print(" +++ using real mode ", real_m)

        real_mh = real_m + 0.5

        # Now things are easy...

        bucket_a_num = np.array(real_mh - 0.5, dtype=int)
        bucket_b_num = np.array(real_mh + 0.5, dtype=int)

        bucket_a_weight = 1.0 - (real_mh - 0.5 - bucket_a_num)
        bucket_b_weight = real_mh + 0.5 - bucket_b_num

        if verbose:
            print(" +++ bucket ", bucket_a_num, " gets ", bucket_a_weight)
            print(" +++ bucket ", bucket_b_num, " gets ", bucket_b_weight)

        assert(np.greater_equal(bucket_a_num, 0).all())
        assert(np.equal(bucket_b_num, bucket_a_num+1).all())

        assert(np.greater_equal(bucket_a_weight, 0).all())
        assert(np.less_equal(bucket_a_weight, 1).all())

        assert np.allclose(bucket_a_weight + bucket_b_weight, 1.0)
        return bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight



    def _modes_to_wavelengths(self, modes):
        return np.pi * 2.0 * self.rsphere / modes


    def _plot_spectrum(self, ke):
        import matplotlib.pyplot as plt

        print(" + resolution in physical space: ", ke.shape)
        ke_spec = np.fft.rfft(np.copy(ke), axis=1)

        print(" + resolution after longitudinal FT transformation: ", ke_spec.shape)

        def _spec_to_buckets_iter(mode_numbers, modes_data, buckets):
            """
            Iterate over all Fourier modes
            m here relates to the number of waves
            """
            for m in range(len(mode_numbers)):
                # Compute real mode (including shortening by being closer to poles)
                # There would be additional number of waves, hence we need to divide by this
                bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = self._ke_spectrum_dist_bucket(mode_numbers[m])

                ampl = np.abs(modes_data[m])

                if bucket_a_num < len(buckets):
                    buckets[bucket_a_num] += ampl

                    if bucket_b_num < len(buckets):
                        buckets[bucket_b_num] += ampl



        def _spec_to_buckets_fast(
                mode_numbers,   # Number of mode. Likely not integer
                modes_data,     # Modal data
                buckets         # Buckets to accumulate
        ):

            bucket_a_num_, bucket_a_weight_, bucket_b_num_, bucket_b_weight_ = self._ke_spectrum_dist_bucket_array(mode_numbers)
            ampl_ = np.abs(modes_data)

            for m in range(len(buckets)):
                if bucket_a_num_[m] < len(buckets):
                    buckets[bucket_a_num_[m]] += ampl_[m]

                    if bucket_b_num_[m] < len(buckets):
                        buckets[bucket_b_num_[m]] += ampl_[m]



        num_buckets = ke_spec.shape[1]-1
        print(" + setting up ", num_buckets, "spectral buckets")
        buckets = np.zeros(num_buckets)

        #
        # Iterate over all longitude stripes
        #
        print(" + processing all latitudes")

        assert(ke_spec.shape[0] == self.sh.lats.shape[0])
        weight_sum = 0
        for i in range(self.sh.lats.shape[0]):
            sys.stdout.write(".")
            sys.stdout.flush()
            #print(" + lat: ", self.sh.lats[i])

            #
            # Compute weight based on width of strip
            #
            # The latitudes are sorted starting with positive ones, then continuing to the negative ones
            #
            if i == 0:
                strip_start = np.pi*0.5
                strip_end = 0.5*(self.sh.lats[0]+self.sh.lats[1])

            elif i == self.sh.lats.shape[0]-1:
                strip_start = 0.5*(self.sh.lats[-1]+self.sh.lats[-2])
                strip_end = -np.pi*0.5

            else:
                strip_start = 0.5*(self.sh.lats[i]+self.sh.lats[i-1])
                strip_end = 0.5*(self.sh.lats[i]+self.sh.lats[i+1])

            weight = (strip_start - strip_end)/np.pi
            #print(" ++ weight: "+str(weight))

            # Compute scalar to multiply modal number with
            # scaling factor \in [
            s = np.cos(self.sh.lats[i])

            m_ = range(num_buckets)
            real_m_ = np.array(m_) / s

            if debug_active:
                a = np.zeros_like(buckets)
                _spec_to_buckets_iter(real_m_, ke_spec[i]*weight, a)

                b = np.zeros_like(buckets)
                _spec_to_buckets_fast(real_m_, ke_spec[i]*weight, b)

                assert np.allclose(a, b)


            _spec_to_buckets_fast(real_m_, ke_spec[i]*weight, buckets)

            weight_sum += weight

        np.allclose(weight_sum, np.pi)


        #
        # Plot results
        #

        import matplotlib.pyplot as plt


        def _plot_buckets(buckets, label=None):
            # Compute wavelengths
            _ = np.arange(len(buckets))
            _[0] = -1.0 # avoid div/0
            wavelengths = self._modes_to_wavelengths(_)

            # bin first and last mode
            if label == None:
                plt.plot(wavelengths[1:-1], buckets[1:-1])
            else:
                plt.plot(wavelengths[1:-1], buckets[1:-1], label=label)

            plt.gca().invert_xaxis()
            plt.xlabel("Wavelength")
            plt.ylabel("Amplitude")

            plt.xscale("log")
            plt.yscale("log")


        #_plot_buckets(buckets, "test")
        _plot_buckets(buckets)


    def _plot_k_lines(self, ax):


        k_ = np.array([100, 200])
        wl_ = self._modes_to_wavelengths(k_)

        #
        # k^-3
        #
        s = 1e15
        y_ = np.power(k_, -3.0)*s
        line = plt.plot(wl_, y_, linestyle="solid", color="gray")
        x = line[0].get_xdata()[1]
        y = line[0].get_ydata()[1]

        ax.annotate(
            "k^-3",
            xy=(x * 1.05, y * 0.35),
            color=line[0].get_color(),
            size=10,
        )

        if False:
            #
            # k^-(5/3)
            #
            s *= np.power(k_[0], -3.0)/np.power(k_[0], -5.0/3.0)
            y_ = np.power(k_, -5.0/3.0)*s
            line = plt.plot(wl_, y_, linestyle="solid", color="gray")
            x = line[0].get_xdata()[1]
            y = line[0].get_ydata()[1]

            ax.annotate(
                "k^-5/3",
                xy=(x * 1.05, y * 1.35),
                color=line[0].get_color(),
                size=10,
            )



    def plot_kinetic_energy_spectrum(self):
        print("Plotting spectrum of kinetic energy")
        """
        Compute
        Ke = 1/2 * m * V^2
        """

        #
        # Compute u and v, prepared for anti-aliasing
        #
        u_phys_data, v_phys_data = self.sh_aa.vrtdiv2uv(self.vrt_spec.data_spectral, self.div_spec.data_spectral)

        #
        # Compute
        #    u*u + v*v
        # and apply anti-aliasing
        #
        V2_phys_data = self.sh_aa.spec2phys(self.sh_aa.phys2spec(u_phys_data * u_phys_data + v_phys_data * v_phys_data))

        #
        # Get mass (which we relate to the height of the SWE)
        # m = phi_pert / g + h0
        #
        m_phys_data = (self.sh_aa.spec2phys(self.phi_pert_spec.data_spectral) / self.grav + self.h0)

        #
        # Finish computation of
        # Ke = 1/2 * m * V^2
        # in physical space
        #
        ke_phys_data = 0.5 * m_phys_data * V2_phys_data

        # Spectral data to get rid of aliasing modes
        ke_spec_data = self.sh_aa.phys2spec(ke_phys_data)

        # Back to physical space
        ke_phys_data = self.sh.spec2phys(ke_spec_data)

        # Pseudo input filename
        input_filename = self.phi_pert_file.replace("phi_pert", "spectrum_kinetic_energy")

        title = input_filename[:]
        title = title.replace('.sweet', '')
        title = title.replace('output_prog_', '')
        _, title = os.path.split(title)

        output_filename = input_filename
        output_filename = output_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')
        output_filename = output_filename.replace('output_prog_', '')

        plt.close()
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        plt.title(title)
        #plt.tight_layout()

        self._plot_spectrum(ke_phys_data)

        self._plot_k_lines(self.ax)

        #plt.legend()
        plt.savefig(output_filename)



    def plot_phi_pert_spectrum(self):
        print("Plotting spectrum of geopotential perturbation")

        # Pseudo input filename
        input_filename = self.phi_pert_file.replace("phi_pert", "spectrum_phi_pert")

        title = input_filename[:]
        title = title.replace('.sweet', '')
        title = title.replace('output_prog_', '')
        _, title = os.path.split(title)

        output_filename = input_filename
        output_filename = output_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')
        output_filename = output_filename.replace('output_prog_', '')


        plt.close()
        plt.title(title)
        plt.tight_layout()
        self.fig, self.ax = plt.subplots(figsize=(6, 4))

        self._plot_spectrum(self.sh.spec2phys(self.phi_pert_spec.data_spectral))

        #plt.legend()
        plt.savefig(output_filename)
 


    def plot_vrt_spectrum(self):
        print("Plotting spectrum of vorticity")

        # Pseudo input filename
        input_filename = self.vrt_file.replace("vrt", "spectrum_vrt")

        title = input_filename[:]
        title = title.replace('.sweet', '')
        title = title.replace('output_prog_', '')
        _, title = os.path.split(title)

        output_filename = input_filename
        output_filename = output_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')
        output_filename = output_filename.replace('output_prog_', '')

        plt.close()
        plt.title(title)
        plt.tight_layout()
        self.fig, self.ax = plt.subplots(figsize=(6, 4))

        self._plot_spectrum(self.sh.spec2phys(self.vrt_spec.data_spectral))

        #plt.legend()
        plt.savefig(output_filename)
 


    def plot_div_spectrum(self):
        print("Plotting spectrum of divergence")

        # Pseudo input filename
        input_filename = self.div_file.replace("div", "spectrum_div")

        title = input_filename[:]
        title = title.replace('.sweet', '')
        title = title.replace('output_prog_', '')
        _, title = os.path.split(title)

        output_filename = input_filename
        output_filename = output_filename.replace('.sweet', '.pdf')
        output_filename = output_filename.replace('/output', '/plot_output')
        output_filename = output_filename.replace('output_prog_', '')

        plt.close()
        plt.title(title)
        plt.tight_layout()
        self.fig, self.ax = plt.subplots(figsize=(6, 4))

        self._plot_spectrum(self.sh.spec2phys(self.div_spec.data_spectral))

        #plt.legend()
        plt.savefig(output_filename)
 
