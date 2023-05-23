import numpy as np
import matplotlib.pyplot as plt
import sys
from mule.SWEETRuntimeParametersScenarios import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from mule.postprocessing.SphereDataSpectral import SphereDataSpectral
import mule.postprocessing.SphereDataOperators as SphereDataOperators

from common import *

import matplotlib.pylab as pylab
params = {'legend.fontsize': 16,
          'legend.title_fontsize': 16,
          'axes.labelsize': 18,
          'axes.titlesize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18}
pylab.rcParams.update(params)

plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": [r'\usepackage{amsmath}']
})


class spectrumSWEET:

    def __init__(self, jobs_data, geometry, fine_key, ref_key):
        self.jobs_data = jobs_data;
        self.geometry = geometry;
        self.fine_key = fine_key;
        self.ref_key = ref_key;
        if self.geometry == "plane":
            self.timescale = 1;
        elif self.geometry == "sphere":
            self.timescale = 1. / (60. * 60.);
        self.rsphere = 6371220

        self.spectra = {}

    def computeKESpectrum(self, t, niters = None):

        """
        See Eq. (5)
        Koshyk, J. N., & Hamilton, K. (2001).
        The horizontal kinetic energy spectrum and spectral budget simulated by a high-resolution troposphere-stratosphere-mesosphere
        GCM. Journal of the Atmospheric Sciences, 58(4), 329–348. https://doi.org/10.1175/1520-0469(2001)058<0329:THKESA>2.0.CO;2
        """
        #self.vrt_spec.data_spectral
        #self.div_spec.data_spectral

        t_str = getFormattedDate(t, self.timescale)

        for job in self.jobs_data.keys():

            self.spectra[job] = {};

            path = getJobPath(self.jobs_data, job)

            if job in [self.fine_key, self.ref_key] or niters is None:
                niters_range = [0]
            else:
                niters_range = niters

            for niter in niters_range:

                if job in [self.fine_key, self.ref_key] or niters is None:
                    vrt_filename = path + "/output_prog_vrt_t" + t_str + ".sweet";
                    div_filename = path + "/output_prog_div_t" + t_str + ".sweet";
                else:
                    vrt_filename = path + "/output_prog_vrt_t" + t_str + "_" + getFormattedIteration(niter) + ".sweet";
                    div_filename = path + "/output_prog_div_t" + t_str + "_" + getFormattedIteration(niter) + ".sweet";


                if not os.path.exists(vrt_filename):
                    spectrum = None
                else:
                    vrt = SphereDataSpectral(vrt_filename, setup_physical=False)
                    div = SphereDataSpectral(div_filename, setup_physical=False)

                    sh = SphereDataOperators.SphereDataOperators(rsphere=self.rsphere)
                    sh.setup_from_file_info(vrt.file_info, anti_aliasing=False)

                    ###def getArrayIndexByModes(n, m):
                    ###    return (m*(2*self.sh.ntrunc-m+1)>>1)+n

                    ###spectrum = np.zeros(self.sh.lats.shape[0])
                    spectrum = np.zeros(sh.lats.shape[0])
                    for m in range(sh.ntrunc):
                        idx = getArrayIndexByModes(m, m, sh.ntrunc)
                        for n in range(m, sh.ntrunc):
                            v = vrt.data_spectral[idx]
                            d = div.data_spectral[idx]

                            if n != 0:
                                spectrum[n] += np.real(1/4*self.rsphere**2/(n*(n+1))*(v*np.conj(v) + d*np.conj(d)))

                            idx += 1

                self.spectra[job][niter] = spectrum

    def modesToWavelengths(self, modes):
        return np.pi * 2.0 * self.rsphere / modes

    ## plot kinetic energy spectra
    def plotKESpectrum(self, title, dirname, output_filename, legend_vars, PinT = False, niters = None, key_fine = None, reference_slopes = True):

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        spectra = self.spectra;

        fontsize=18
        figsize=(9, 7)

        plt.figure(1, figsize=figsize)

        ax = plt.gca()

        if reference_slopes:

            ###ymin, ymax = np.log([spec[1:].min(), spec[1:].max()])
            ymin, ymax = np.log([1e0,1e4])
            xmin, xmax = np.log([2e2, 6e2])

            for slope in (-5./3., -3.):
                y1 = ymin
                y2 = slope * (xmax - xmin) + ymin
                ax.loglog(np.exp([xmax, xmin]), np.exp([y1, y2]), 'black', linestyle = '--')
                if slope == -3.:
                    ax.annotate(r'$n^{-3}$', xy = (np.exp(xmin), np.exp(y2 - 2)), fontsize = fontsize)
                else:
                    ax.annotate(r'$n^{-5/3}$', xy = (np.exp(xmin), np.exp(y2 + 1)), fontsize = fontsize)

            ######## TODO
            ######## TODO: Plot reference slopes
            ######## TODO
            #######en_ref53=np.array([])
            #######iref=(spectra[self.ref_key][0][1]/50.0)/np.power(float(1), -float(5.0/3.0))	
            #######for tmp in r_ref3:
            #######    ytmp=np.power(tmp, -float(5.0/3.0))*iref
            #######    en_ref53=np.append(en_ref53, [ytmp])

            #######    r_ref3_len=xL_max*1000/r_ref3[1:]
            #######    #plt.loglog(r_ref53, en_ref53, '-.', color='black')
            #######    plt.loglog(r_ref3_len, en_ref3[1:], '-.', color='black')
            #######    plt.loglog(r_ref3_len, en_ref53[1:], '-.', color='black')

            #######    ax.annotate("$k^{-5/3}$", xy=(r_ref3_len[-1]-10, en_ref53[-1]), fontsize=fontsize)
            #######    ax.annotate("$k^{-3}$", xy=(r_ref3_len[-1]-10, en_ref3[-1]), fontsize=fontsize)

            #######    ax.spines['top'].set_visible(False)
            #######    ax.spines['right'].set_visible(False)
            #######    ax.xaxis.set_label_coords(0.5, -0.075)
            #######    ax.set_facecolor('xkcd:white')

        # invert axis for wavelength
        ax.invert_xaxis()

        ax.xaxis.set_label_coords(0.5, -0.075)

        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)

        #plt.xticks(labelsx, fontsize=fontsize)
        plt.xlabel("Horizontal wavelength ($km$)", fontsize=fontsize)

        #plt.yticks(labelsy, fontsize=fontsize)
        plt.ylabel("Kinetic Energy ($m^2 s^{-2}$)", fontsize=fontsize)

        markers = ['']
        ##markers = ['+', 'o', 'x', 'D', 's', '^', 'v', '1', '2', '3', '4'];
        linestyles = ['-', '--', ':', '-.']
        ##linestyles = ['--', ':', '-.']

        ## sort jobs
        list_jobs = list(spectra.keys())
        plot_fine = False
        plot_ref = False
        if self.ref_key in list_jobs:
            list_jobs.remove(self.ref_key)
            plot_ref = True
        if self.fine_key in list_jobs:
            list_jobs.remove(self.fine_key)
            plot_fine = True

        jobs_to_plot = sortJobs(list_jobs, legend_vars, self.jobs_data)

        if plot_ref:
            jobs_to_plot.append(self.ref_key)
        if plot_fine:
            jobs_to_plot.append(self.fine_key)

        if not PinT:
            c = 0
            for job in jobs_to_plot:
                if not spectra[job] is None:
                    ###y = spectra[job]["spectrum"];
                    y = spectra[job][0];
                    if y is None:
                        continue
                    x = np.arange(len(y))
                    x[0] = -1.0
                    x = self.modesToWavelengths(x) / 1e3
                    l = getLabel(legend_vars, job, self.jobs_data);
                    markevery = int(len(x) / 10 + 2 * c);

                    plt.loglog(x[1:], y[1:], markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)], markevery = markevery, markersize = 5, label=l)	
                    c = c + 1

        else:
            ijob = 0;
            for job in jobs_to_plot:
                initer = 0;
                ####if job == self.ref_key and self.fine_key == self.ref_key:
                ####    continue
                if job in [self.fine_key, self.ref_key]:
                    niters_range = [0];
                else:
                    niters_range = niters;
                for niter in niters_range:
                    y = spectra[job][niter];
                    if y is None:
                        continue
                    x = np.arange(len(y))
                    x[0] = -1.0
                    x = self.modesToWavelengths(x) / 1e3
                    markevery = int(len(y) / 10);
                    markersize = 4.
                    if job == self.ref_key:
                        l = "ref";
                    elif job == self.fine_key:
                        l = "fine";
                    else:
                        l = getLabel(legend_vars, job, self.jobs_data) + "; " + r'$k = {}$'.format(niter);

                    markers = ["", 's', 'o', 'D']

                    color = getPlotAttribute("color", ijob)
                    marker = markers[initer % len(markers)]
                    linestyle = getPlotAttribute("linestyle", initer)

                    ###plt.loglog(x[1:], y[1:], markers[c % len(markers)], linestyle=linestyles[c % len(linestyles)],  markevery = markevery, label=l)	
                    plt.loglog(x[1:], y[1:], color = color, marker = marker, linestyle = linestyle,  markevery = markevery, markersize = markersize, label = l)
                    initer += 1

                ijob += 1

        if title != '':
            plt.title(title, fontsize=fontsize)

        plt.legend(title = setTitleLegend(legend_vars))

        plt.savefig(dirname + "/" + output_filename + ".pdf", transparent=True, bbox_inches='tight', pad_inches=0)

        plt.close()
