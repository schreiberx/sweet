import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d, CubicSpline
from common import *
import seaborn as sns
from matplotlib.colors import LogNorm, Normalize

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


class errorSWEET:

    def __init__(self, jobs_data, jobs_solution, ref_key, ref_type, geometry, var_error = "prog_phi_pert"):

        self.jobs_data = jobs_data;
        self.jobs_solution = jobs_solution;
        self.ref_key = ref_key;
        self.ref_type = ref_type
        self.geometry = geometry;
        self.var_error = var_error;

        if not (ref_type == "ref" or ref_type == "fine"):
            sys.exit("Wrong ref_type: " + ref_type);

        if self.geometry == "plane" or self.geometry == "scalar":
            self.timescale = 1;
        elif self.geometry == "sphere":
            self.timescale = 1. / (60. * 60.);

        self.keys_errors = ["err_L1", "err_L2", "err_Linf"];

        self.err = {};
        self.err_time_iters_jobs = {};
        self.err_iters_time = {};
        self.err_time_iters = {};
        self.err_iters = {};
        self.err_time_iters_array = {};
        self.err_iters_array = {};
        self.err_iters_threshold = {}


        ## Residuals
        self.residuals = {};


    def readPinTErrorsFromFile(self, job, t, niter, spectral = False, pint_type = "xbraid"):

        fdate = t;

        ## get reference data (to compute relative error)
        if self.geometry == "scalar":
            ref_L1 = 1;
            ref_L2 = 1;
            ref_Linf = 1;
        else:
            data_ref = self.jobs_solution[self.ref_key][fdate];
            if data_ref is None:
                return None;
            shape_ref = data_ref.shape;
            if len(shape_ref) == 1:
                nx = 1;
                ny = shape_ref[0];
            else:
                (ny, nx) = shape_ref;
            ####ref_L1 = np.sum(np.absolute(data_ref))/(nx * ny);
            ####ref_L2 = math.sqrt( np.sum(np.square(data_ref))  / (nx * ny)  );
            ####ref_Linf = np.max(np.absolute(data_ref))
            ref_L1 = 1.;
            ref_L2 = 1.;
            ref_Linf = 1.;

        err = {};

        ## physical or spectral solution
        spec_str = "";
        if spectral:
            spec_str = "spec_";

        ## errors computed w.r.t. ref or fine solution
        ref_str = self.ref_type;

        ## check if file exists
        path = os.path.basename(self.jobs_data[job]["jobgeneration.p_job_dirpath"]);
        fname = path + "/" + pint_type + "_error_" + spec_str + ref_str + "_" + self.var_error + "_t" + t + "_" + getFormattedIteration(niter) + ".csv";
        if not os.path.isfile(fname):
            return None;

        ## read errors
        lines = [line.rstrip() for line in open(fname)];

        spl = lines[0].split();
        assert spl[0] == "#BASESOLUTION";
        assert spl[1] == ref_str;

        spl = lines[1].split();
        assert spl[0] == "#VAR";
        assert spl[1] == self.var_error, (spl[1], self.var_error);

        spl = lines[2].split();
        assert spl[0] == "#ITERATION";
        assert int(spl[1]) == niter;

        spl = lines[3].split();
        assert spl[0] == "#TIMESLICE";

        spl = lines[4].split();
        assert spl[0] == "#TIMEFRAMEEND";
        assert np.abs(float(spl[1]) - float(t)) < 1e-7, (spl[1], t, fname);

        ## read L1, L2, Linf errors
        if not spectral:
            if len(lines) >= 9:
                for i in range(5, len(lines)):
                    spl = lines[i].split();
                    err[spl[0]] = float(spl[1]);
            else:
                spl = lines[5].split();
                assert spl[0] == "errL1";
                err["err_L1"] = float(spl[1]) / ref_L1;

                spl = lines[6].split();
                assert spl[0] == "errL2";
                err["err_L2"] = float(spl[1]) / ref_L2;

                spl = lines[7].split();
                assert spl[0] == "errLinf";
                err["err_Linf"] = float(spl[1]) / ref_Linf;
        ## read Linf(Rnorm) errors
        else:
            self.keys_errors = [];
            for iline in range(5, len(lines)):
                spl = lines[iline].split();
                assert len(spl) == 3;
                assert spl[0] == "errLinf";
                err["rnorm_" + spl[1]] = float(spl[2]);
                if int(spl[1]) <= 8:
                    continue;
                self.keys_errors.append("rnorm_" + spl[1]);

        return err;


    def readPinTResiduals(self, pint_type = "xbraid"):

        for job in self.jobs_data.keys():
            if job == self.ref_key or job == "ref":
                continue;
            if not self.jobs_data[job]["runtime." + pint_type + "_enabled"]:
                continue;
            self.residuals[job] = np.empty((0, 2), dtype = float);
            niter = 1;
            while True:
                r = self.readPinTResidualsJobNiter(job, niter);
                if r == None:
                    break;
                self.residuals[job] = np.vstack((self.residuals[job], np.array([niter, r])));
                niter += 1;


    ## read residual at given iteration
    def readPinTResidualsJobNiter(self, job, niter):

        fname = getJobPath(self.jobs_data, job) + "/residual" + "_" + getFormattedIteration(niter) + ".csv";

        if not os.path.isfile(fname):
            return None;

        lines = [line.rstrip() for line in open(fname)];

        spl = lines[0].split();
        assert spl[0] == "#SWEET";

        spl = lines[1].split();
        assert spl[0] == "#FORMAT";
        assert spl[1] == "ASCII";

        spl = lines[2].split();
        assert spl[0] == "#PRIMITIVE";
        assert spl[1] == "SCALAR";

        spl = lines[3].split();
        assert len(spl) == 1;

        if math.isnan(float(spl[0])):
            return None;

        return float(spl[0]);


    def computeErrorsPinT(self, unify_time_steps = False, read_errors_from_file = False, pint_type = "parareal", spectral = False):

        ## get instants in which the error is computed
        tts = list(self.jobs_solution[self.ref_key].keys());

        ## get errors for each time and iteration --> dict[t][niter][job]
        for t in tts:
            niter = 0;
            self.err_time_iters_jobs[t] = {};
            while True:
                if not read_errors_from_file:
                    sys.exit("Offline error computation not available");
                else:
                    tmp = {};
                    found_simulation = False;
                    self.err_time_iters_jobs[t][niter] = {};
                    for job in self.jobs_solution.keys():
                        ## no error in reference simulation
                        if job == self.ref_key or job == "ref":
                            continue;
                        ## no error if not PinT simulation
                        if not self.jobs_data[job]["runtime." + pint_type + "_enabled"]:  ## not parareal
                            continue;
                        e = self.readPinTErrorsFromFile(job, t, niter, spectral, pint_type);
                        if not e is None:
                            self.err_time_iters_jobs[t][niter][job] = e;
                            found_simulation = True;
                    niter += 1;
                    ## no more jobs available
                    if not found_simulation:
                        break;

        ## reorder err_time_iters --> dict[job][t][niter] and dict[job][niter][t]
        tmp = {};
        for t in tts:
            for niter in self.err_time_iters_jobs[t].keys():
                for job in self.jobs_solution.keys():
                    if job == self.ref_key or job == "ref":
                        continue;
                    if not self.jobs_data[job]["runtime." + pint_type + "_enabled"]:
                        continue;
                    if job in self.err_time_iters_jobs[t][niter].keys():
                        if not job in self.err_time_iters.keys():
                            self.err_time_iters[job] = {};
                        if not t in self.err_time_iters[job].keys():
                            self.err_time_iters[job][t] = {};
                        if not job in self.err_iters_time.keys():
                            self.err_iters_time[job] = {};
                        if not niter in self.err_iters_time[job].keys():
                            self.err_iters_time[job][niter] = {};
                        self.err_time_iters[job][t][niter] = self.err_time_iters_jobs[t][niter][job];
                        self.err_iters_time[job][niter][t] = self.err_time_iters_jobs[t][niter][job];

        ## get maximum error in time per iteration
        for job in self.err_iters_time.keys():
            self.err_iters[job] = {};
            for niter in self.err_iters_time[job].keys():
                max_error = {};
                self.err_iters[job][niter] = {};
                for key in self.keys_errors:
                    self.err_iters[job][niter][key] = -1e10;
                    for t in self.err_iters_time[job][niter].keys():
                        ## error = 0 at t = 0
                        if float(t) == 0:
                            continue;
                        ## update maximum error
                        if self.err_iters_time[job][niter][t][key] > self.err_iters[job][niter][key]:
                            self.err_iters[job][niter][key] = self.err_iters_time[job][niter][t][key];
                        ## if diverged -> stop
                        if float(t) > 0 and (np.isnan(self.err_iters_time[job][niter][t][key]) or np.isinf(self.err_iters_time[job][niter][t][key])):
                            self.err_iters[job][niter][key] = np.nan;
                            break;

        ## convert dicts to array (for plotting)
        self.convertErrorDictToArray();

        ## get number of iterations for reaching a given threshold
        self.getIterationsForErrorThreshold();

    def convertErrorDictToArray(self):

        ## arrays with max error in time per iteration; array[:, 0] = niter; array[:, 1] = error(niter)
        self.err_iters_array = {};
        for job in self.err_iters.keys():
            self.err_iters_array[job] = {};
            for key in self.keys_errors:
                self.err_iters_array[job][key] = np.empty((0, 2), dtype = float);
                for niter in self.err_iters[job].keys():
                    self.err_iters_array[job][key] = np.vstack((self.err_iters_array[job][key], np.array([niter, self.err_iters[job][niter][key]])));

        ## arrays wit error per time; array[:, 0] = time; array[:, 1] = error(time)
        self.err_time_iters_array = {};
        for job in self.err_time_iters.keys():
            self.err_time_iters_array[job] = {};
            for key in self.keys_errors:
                self.err_time_iters_array[job][key] = {};
                for niter in self.err_time_iters[job].keys():
                    self.err_time_iters_array[job][key][niter] = np.empty((0, 2), dtype = float);
                    for t in self.err_time_iters[job][niter].keys():
                        self.err_time_iters_array[job][key][niter] = np.vstack((self.err_time_iters_array[job][key][niter], np.array([float(t), self.err_time_iters[job][niter][t][key]])));

    def getIterationsForErrorThreshold(self, read_from_file = False, pint_type = "xbraid"):

        ## information already computed during SWEET execution
        if read_from_file:
            for job in self.jobs_data.keys():
                self.readIterationsForErrorThresholdsFromFile(job, pint_type);
            return;

        ## else
        thresholds = np.arange(0, -17, -1, dtype = float);
        for job in self.err_iters.keys():
            self.err_iters_threshold[job] = {};
            for key in self.err_iters_array[job].keys():
                self.err_iters_threshold[job][key] = {};
                for threshold in thresholds:
                    eps = 10 ** threshold;
                    self.err_iters_threshold[job][key][eps] = 1e16;
                    for niter in range(self.err_iters_array[job][key].shape[0]):
                        if self.err_iters_array[job][key][niter, 1] < eps and niter < self.err_iters_threshold[job][key][eps]:
                            self.err_iters_threshold[job][key][eps] = niter;

    ## read number of iterations for reaching error thresholds from output files
    def readIterationsForErrorThresholdsFromFile(self, job, pint_type = "xbraid"):

        ref_str = self.ref_type;

        path = getJobPath(self.jobs_data, job)
        fname = path + "/" + pint_type + "_iter_thresholds_" + ref_str + "_" + self.var_error + ".csv";

        if not os.path.isfile(fname):
            return None;

        lines = [line.rstrip() for line in open(fname)];

        spl = lines[0].split();
        assert spl[0] == "#BASESOLUTION";
        assert spl[1] == ref_str;

        spl = lines[1].split();
        assert spl[0] == "#VAR";
        assert spl[1] == self.var_error, (spl[1], self.var_error);

        self.err_iters_threshold[job] = {};
        for iline in range(2, len(lines)):
            spl = lines[iline].split();
            if (len(spl) < 3):
                continue;
            key = spl[0];
            eps = np.power(10, -float(spl[1]));
            niter = int(spl[2]);
            if key not in self.err_iters_threshold[job].keys():
                self.err_iters_threshold[job][key] = {};
            self.err_iters_threshold[job][key][eps] = niter;

    def removeJobs(self, delete_vars, delete_vals):

        jobs_to_remove = []

        ivar = 0
        for var in delete_vars:
            for job in self.err_iters_array.keys():
                if job in jobs_to_remove:
                    continue
                if self.jobs_data[job][var] in delete_vals[ivar]:
                    jobs_to_remove.append(job)
                    continue
            ivar += 1

        for job in jobs_to_remove:
            del(self.err_iters_array[job])
            if job in self.residuals.keys():
                del(self.residuals[job])


    ## plot error or residual along iterations for a given set of simulation and a given error norm
    def plotErrorPinTAlongIterations(self, dirname, plot_type, err_type, legend_vars, groups_vars, common_plot_attributes, filter_vars = None, filter_vals = None, relative_error = True, max_iter = None, plot_legend = True, ncol_legend = 1, first_plot_idx = 0, ylim = None, linewidth = .75, legend_fontsize = 16):

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        if plot_type == "error":
            array_jobs = self.err_iters_array.copy();
        elif plot_type == "residual":
            array_jobs = self.residuals.copy();
        else:
            sys.exit("Wrong plot_type: " + plot_type);

        if self.ref_key in array_jobs:
            array_jobs.remove(self.ref_key)

        jobs_to_plot = sortJobs(array_jobs, legend_vars, self.jobs_data)

        jobs_to_plot = filterJobs(jobs_to_plot, filter_vars, filter_vals, self.jobs_data)

        ## group jobs based on given parameters
        groups = createGroups(jobs_to_plot, groups_vars, self.jobs_data);

        fig, ax = plt.subplots();
        igroup = first_plot_idx;
        for group in groups.values():

            ijob = 0;
            for job in group:

                if plot_type == "error":
                    array = array_jobs[job][err_type];
                elif plot_type == "residual":
                    array = array_jobs[job];
                if (array.shape[0] == 0):
                    continue;

                ## remove too large errors (instabilities)
                array2 = np.empty((0,2), dtype = float)
                for i in range(array.shape[0]):
                    if array[i, 1] < 1e10:
                        array2 = np.vstack((array2, array[i, :]))
                    else:
                        break
                array = array2

                ## get color, marker, linestyle
                cml = {};
                for pa in plot_attributes.keys():
                    if pa in common_plot_attributes:
                        cml[pa] = getPlotAttribute(pa, igroup);
                    else:
                        cml[pa] = getPlotAttribute(pa, ijob);

                if max_iter is None:
                    range_iter = array.shape[0];
                else:
                    range_iter = min(max_iter + 1, array.shape[0])

                label = getLabel(legend_vars, job, self.jobs_data);
                ax.plot(array[:range_iter, 0], array[:range_iter, 1], label = label, color = cml["color"], marker = cml["marker"], linestyle = cml["linestyle"], linewidth = linewidth);
                ijob += 1;

            igroup += 1;

        if not (ylim is None):
            ax.set_ylim(ylim)

        ax.set_yscale("log");
        ax.set_xlabel("Iteration");
        if relative_error:
            if "rnorm" in err_type:
                ax.set_ylabel("Relative error " + getPrettyName(err_type + "_" + self.var_error) );
            else:
                ax.set_ylabel("Relative error " + getPrettyName(err_type) );
        else:
            ax.set_ylabel("Error " + getPrettyName(err_type) );
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        if plot_legend:
            ax.legend(title = setTitleLegend(legend_vars), ncol = ncol_legend, loc='best', fontsize = legend_fontsize);

        plt.tight_layout();

        figname = dirname + "/" + plot_type + "_" + self.ref_type + "_" + err_type;
        if not filter_vars is None:
            figname = setFilenameFilterVars(figname, filter_vars, filter_vals);

        plt.savefig(figname + ".pdf");
        plt.close();


    ## compute difference between reference and (one or more) solutions at given time
    def computeDiff(self, i_data, t, spectral = False):

        fdate = getFormattedDate(t, self.timescale);

        data_ref = i_data[self.ref_key][fdate];

        ## find the lowest spatial resolution among all jobs
        first = True;
        nx_first = 0;
        ny_first = 0;
        for job in i_data.keys():
            ## reference solution: nothing to do
            if job == self.ref_key:
                continue;
            ## solution is not available at this time (maybe due to instabilities): nothing to do
            if i_data[job][fdate] is None:
                continue;
            ## first found solution
            if first:
                ## get spatial dimensions
                shape = i_data[job][fdate].shape;
                if len(shape) == 1:
                    nx_first = 1;
                    ny_first = shape[0];
                else:
                    (ny_first, nx_first) = shape;
                job_low_resolution = job;
                first = False;
            else:
                shape = i_data[job][fdate].shape;
                if len(shape) == 1:
                    nx_cmp = 1;
                    ny_cmp = shape[0];
                else:
                    (ny_cmp, nx_cmp) = shape;

                if nx_cmp < nx_first and ny_cmp < ny_first:
                    nx_first = nx_cmp;
                    ny_first = ny_cmp;
                    job_low_resolution = job;

                if ny_cmp == 0:
                    print (nx_cmp, ny_cmp, i_data[job][fdate], job, self.data["path"][job]);

        nx_cmp = nx_first;
        ny_cmp = ny_first;

        ## interpolate ref to the lowest resolution
        if not spectral:
            data_ref_low = self.interpolate(data_ref, ny_cmp, nx_cmp);
        else:
            data_ref_low = self.spectralTruncate(data_ref, ny_cmp, nx_cmp);

        ## interpolate jobs to the lowest resolution and compute difference w.r.t. the reference solution
        data_diff = {};
        for job in i_data.keys():
            if job == self.ref_key:
                continue;
            if i_data[job][fdate] is None:
                data_diff[job] = None;
                continue;

            if not spectral:
                data_low_res = self.interpolate(i_data[job][fdate], ny_cmp, nx_cmp);
            else:
                data_low_res = self.spectralTruncate(i_data[job][fdate], ny_cmp, nx_cmp);
            data_diff[job] = data_low_res - data_ref_low;

        return data_diff, data_ref_low, nx_cmp, ny_cmp;




    ## compute errors in physical space
    ## interpolate solutions to the lowest spatial resolution and compute difference
    def computePhysicalErrors(self, t, data = None, relative_error = True):

        ## compute errors for all previously read jobs
        if data is None:
            data = self.jobs_solution;
        data_diff, data_ref, nx_cmp, ny_cmp = self.computeDiff(data, t, spectral = False);

        if relative_error:
            shape_ref = data_ref.shape;
            if len(shape_ref) == 1:
                nx_ref = 1;
                ny_ref = shape[0];
            else:
                (ny_ref, nx_ref) = shape_ref;
            ref_L1 = np.sum(np.absolute(data_ref))/(nx_ref * ny_ref);
            ref_L2 = math.sqrt( np.sum(np.square(data_ref))  / (nx_ref * ny_ref)  );
            ref_Linf = np.max(np.absolute(data_ref))
        else:
            ref_L1 = 1;
            ref_L2 = 1;
            ref_Linf = 1;

        for job in data_diff.keys():
            if job == self.ref_key:
                continue;
            if not (job in self.err.keys()):
                self.err[job] = {};

            if data_diff[job] is None:
                self.err[job]["err_L1"] = None;
                self.err[job]["err_L2"] = None;
                self.err[job]["err_Linf"] = None;
                continue;

            self.err[job]["err_L1"] = np.sum(np.absolute(data_diff[job]))/(nx_cmp * ny_cmp) / ref_L1;
            self.err[job]["err_L2"] = math.sqrt( np.sum(np.square(data_diff[job]))  / (nx_cmp * ny_cmp)  ) / ref_L2;
            self.err[job]["err_Linf"] = np.max(np.absolute(data_diff[job])) / ref_Linf

    ## compute errors in spectral space
    ## truncate spectra to the lowest spectral resolution and compute max difference
    def computeSpectralErrors(self, t, data = None, relative_error = True):

        ## compute errors for all previously read jobs
        if data is None:
            data = self.jobs_solution;
        data_diff, data_ref, nx_cmp, ny_cmp = self.computeDiff(data, t, spectral = True);

        for job in data_diff.keys():
            if job == self.ref_key:
                continue;
            if not (job in self.err.keys()):
                self.err[job] = {};

            for rnorm in [16, 32, 64, 128, 256, 512, 1024]:
                self.err[job]["rnorm_" + str(rnorm)] = getMaxAbsRnorm(data_diff[job], rnorm, -1, verbose = False)
                if relative_error and not(data_diff[job] is None):
                    self.err[job]["rnorm_" + str(rnorm)] /= getMaxAbsRnorm(data_ref, rnorm, -1, verbose = False)

    def spectralTruncate(self, data_high_resolution, ny_cmp, nx_cmp):

        ## get high spatial resolution
        shape = data_high_resolution.shape;
        if len(shape) == 1:
            nx_ref = 1;
            ny_ref = shape[0];
        else:
            (ny_ref, nx_ref) = shape;

        data_low_resolution = np.empty((ny_cmp, nx_cmp))

        for m in range(ny_cmp):
            for n in range(m, nx_cmp):
               data_low_resolution[m, n] = data_high_resolution[m, n]

        return data_low_resolution


    ## interpolate solution to the lower spatial resolution (nx_cmp, ny_cmp)
    def interpolate(self, data_high_resolution, ny_cmp, nx_cmp):

        ## get high spatial resolution
        shape = data_high_resolution.shape;
        if len(shape) == 1:
            nx_ref = 1;
            ny_ref = shape[0];
        else:
            (ny_ref, nx_ref) = shape;

        multiplier_j = (ny_ref)/(ny_cmp)
        multiplier_i = (nx_ref)/(nx_cmp)

        if multiplier_i == 1 and multiplier_j == 1: #Grids are the same
            return data_high_resolution;
        elif multiplier_i > 1 and multiplier_j > 1 :
            #Comparison via interpolation
            #print("Interpolation")
            # A-grid REFERENCE (file1) - sweet outputs only A grids physical space
            dx_ref=1.0/(nx_ref)
            dy_ref=1.0/(ny_ref)

            x_ref = np.arange(0, 1, dx_ref)
            x_ref = np.linspace(0, 1, nx_ref, endpoint=False)

            y_ref = np.arange(0, 1, dy_ref)
            y_ref = np.linspace(0, 1, ny_ref, endpoint=False)

            x_ref += dx_ref/2
            y_ref += dy_ref/2
            X_ref, Y_ref = np.meshgrid(x_ref, y_ref)

            #Create cubic interpolation of reference file
            interp_spline = RectBivariateSpline(y_ref, x_ref, data_high_resolution)

            #A-grid cmp file (file2)
            dx_cmp=1.0/nx_cmp
            dy_cmp=1.0/ny_cmp

            x_cmp = np.arange(0, 1, dx_cmp)
            x_cmp = np.linspace(0, 1, nx_cmp, endpoint=False)

            y_cmp = np.arange(0, 1, dy_cmp)
            y_cmp = np.linspace(0, 1, ny_cmp, endpoint=False)

            x_cmp += dx_cmp/2
            y_cmp += dy_cmp/2
            X_cmp, Y_cmp = np.meshgrid(x_cmp, y_cmp)

            #Get reduced reference resolution
            data_low_resolution = interp_spline(y_cmp, x_cmp)

            return data_low_resolution;

        elif multiplier_i >= 1 or multiplier_j >= 1:  ## 1D interpolation

            dy_ref=1.0/(ny_ref)

            y_ref = np.linspace(0, 1, ny_ref, endpoint=False)
            y_ref += dy_ref/2

            #Create cubic interpolation of reference file
            interp_spline = CubicSpline(y_ref, data_high_resolution);

            dy_cmp=1.0/ny_cmp

            y_cmp = np.arange(0, 1, dy_cmp)
            y_cmp = np.linspace(0, 1, ny_cmp, endpoint=False)

            data_low_resolution = interp_spline(y_cmp)

            return data_low_resolution;


    def plotErrorSerial(self, dirname, err_type, x_var, legend_vars, groups_vars, filter_vars = None, filter_vals = None, x_logscale = True, legend_vars_plot_attributes = None):

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        list_jobs = list(self.jobs_data.keys());
        if self.ref_key in list_jobs:
            list_jobs.remove(self.ref_key);

        jobs_to_plot = sortJobs(list_jobs, legend_vars, self.jobs_data)

        jobs_to_plot = filterJobs(jobs_to_plot, filter_vars, filter_vals, self.jobs_data)

        groups = createGroups(jobs_to_plot, groups_vars, self.jobs_data);

        if not (legend_vars_plot_attributes is None):
            common_plot_attributes = {};


        if type(err_type) == str:
            err_types = [err_type]
            several_err_types = False
        else:
            err_types = err_type
            several_err_types = True


        fig, ax = plt.subplots();
        igroup = 0;
        iplot = {"color": 0, "marker": 0,  "linestyle": 0}

        for err_type in err_types:
            for group in groups.values():

                x_plot = np.empty(0, dtype = float);
                y_plot = np.empty(0, dtype = float);

                for job in group:

                    x = self.jobs_data[job][x_var]
                    y = self.err[job][err_type]

                    if y is None:
                        continue;

                    x_plot = np.append(x_plot, x);
                    y_plot = np.append(y_plot, y);

                ## sort x values
                idx = np.argsort(x_plot);
                x_plot = x_plot[idx];
                y_plot = y_plot[idx];

                ## get color, marker, linestyle
                ## each simulation has a different color, marker, linestyle
                cml = {};
                for pa in plot_attributes.keys():
                    cml[pa] = getPlotAttribute(pa, igroup);

                ## jobs with common parameters share plot attributes
                ## if legend_vars = [var1, var2], then e.g.
                ## legend_vars_plot_attributes = [["color"], ["marker", "linestyle"]]
                ## then all jobs with var1 share same color
                ## and all jobs with var2 share same marker and linestyle
                if not (legend_vars_plot_attributes is None or legend_vars is None):
                    ## get legend vals for this job
                    ivar = 0;
                    for var in legend_vars:
                        val = self.jobs_data[job][var]
                        pas = legend_vars_plot_attributes[ivar];  ## array, e.g. ["color"]
                        ## check if a plot attribute has already been set to this val
                        key = var + str(val);
                        if key in common_plot_attributes.keys():
                            for pa in pas:
                                cml[pa] = common_plot_attributes[key][pa]
                        else:
                            common_plot_attributes[key] = {};
                            for pa in pas:
                                cml[pa] = getPlotAttribute(pa, iplot[pa])
                                common_plot_attributes[key][pa] = cml[pa]
                                iplot[pa] += 1
                        ivar += 1

                if several_err_types:
                    label = getPrettyName(err_type + "_" + self.var_error)
                else:
                    label = getLabel(legend_vars, job, self.jobs_data);
                ax.plot(x_plot, y_plot, label = label, color = cml["color"], marker = cml["marker"], linestyle = cml["linestyle"], linewidth = .75);

                igroup += 1;

        if several_err_types or not (legend_vars is None):
            ax.legend();

        if several_err_types:
            ax.set_ylabel("Relative error")
        else:
            title = "Error " + getPrettyName(err_type + "_" + self.var_error);
            ax.set_title(title);
            ax.set_ylabel("Relative error " + getPrettyName(err_type + "_" + self.var_error))

        ax.set_xlabel(getPrettyName(x_var));

        ax.set_yscale("log")
        if (x_logscale):
            ax.set_xscale("log");

        if several_err_types:
            figname = dirname + "/all_err";
        else:
            figname = dirname + "/" + err_type;
        if not (groups_vars is None):
            for var in groups_vars:
                figname += "__" + var;
        figname += ".pdf";

        plt.tight_layout()

        fig.savefig(figname);
        plt.close();



    def plotMapErrorsAtIteration(self, dirname, var_x, var_y, niter, plot_type, err_type, filter_vars = None, filter_vals = None, relative_error = True, axis_labels_power = True, dummy_x_values = [], dummy_y_values = [], vmin = None, vmax = None):

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        if plot_type == "error":
            array_jobs = self.err_iters_array;
        elif plot_type == "residual":
            array_jobs = self.residuals;
        else:
            sys.exit("Wrong plot_type: " + plot_type);

        if self.ref_key in array_jobs:
            array_jobs.remove(self.ref_key)

        jobs_to_plot = filterJobs(list(array_jobs.keys()), filter_vars, filter_vals, self.jobs_data)

        x_array = np.empty(0, dtype = float)
        y_array = np.empty(0, dtype = float)
        xy_dict = {};

        ## group jobs based on given parameters
        fig, ax = plt.subplots();
        for job in jobs_to_plot:

            if plot_type == "error":
                array = array_jobs[job][err_type];
            elif plot_type == "residual":
                array = array_jobs[job];
            if (array.shape[0] == 0):
                continue;

            ## remove too large errors (instabilities)
            array2 = np.empty((0,2), dtype = float)
            for i in range(array.shape[0]):
                if array[i, 1] < 1e10:
                    array2 = np.vstack((array2, array[i, :]))
                else:
                    break
            array = array2

            ## get x, y and z for this job
            x = self.jobs_data[job][var_x]
            y = self.jobs_data[job][var_y]
            x_array = np.append(x_array, x)
            y_array = np.append(y_array, y)
            if niter < array.shape[0]:
                xy_dict[(x, y)] = array[niter, 1]
            else:
                xy_dict[(x, y)] = -1

        ## get min and max z values
        vals = list(xy_dict.values())
        if vmin is None:
            if np.any(np.array(vals) >= 0):
                vmin = np.min([z for z in vals if z >= 0])
            else:
                vmin = 0.
        if vmax is None:
            vmax = np.max(vals)

        x_array_with_dummy = x_array.copy()
        y_array_with_dummy = y_array.copy()
        ## add dummy x and y values
        assert np.all(np.in1d(x_array, dummy_x_values)) == False
        assert np.all(np.in1d(y_array, dummy_y_values)) == False
        for x in dummy_x_values:
            for y in np.append(y_array_with_dummy, dummy_y_values):
                x_array_with_dummy = np.append(x_array_with_dummy, x)
                y_array_with_dummy = np.append(y_array_with_dummy, y)
                ###xy_dict[(x, y)] = 1e20
        for y in dummy_y_values:
            for x in np.append(x_array_with_dummy, dummy_x_values):
                x_array_with_dummy = np.append(x_array_with_dummy, x)
                y_array_with_dummy = np.append(y_array_with_dummy, y)
                ####xy_dict[(x, y)] = 1e20

        ## remove duplicates in x and y
        x_array = np.unique(x_array)
        y_array = np.unique(y_array)
        x_array_with_dummy = np.unique(x_array_with_dummy)
        y_array_with_dummy = np.unique(y_array_with_dummy)

        ## sort x and y
        x_array = np.sort(x_array)
        y_array = np.sort(y_array)
        x_array_with_dummy = np.sort(x_array_with_dummy)
        y_array_with_dummy = np.sort(y_array_with_dummy)

        ## construct z matrix
        z_matrix = np.ones((x_array.size, y_array.size), dtype = float) * (-1)
        for key in xy_dict.keys():
            x = key[0]
            y = key[1]
            z = xy_dict[key]
            ix = np.where(x_array == x)[0][0]
            iy = np.where(y_array == y)[0][0]
            z_matrix[ix, iy] = z

        if axis_labels_power:
            xticklabels = [r'$10^{{{}}}$'.format(int(np.log10(x+1))) for x in x_array]
            yticklabels = [r'$10^{{{}}}$'.format(int(np.log10(y+1))) for y in y_array]
            if 0 in x_array:
                xticklabels[np.argwhere(x_array == 0)[0][0]] = r'$0$'
            if 0 in y_array:
                yticklabels[np.argwhere(y_array == 0)[0][0]] = r'$0$'
        else:
            xticklabels = x_array
            yticklabels = y_array

        ######## set dummy tick to dummy axis values
        ######for x in dummy_x_values:
        ######    xticklabels[np.argwhere(x_array == x)[0][0]] = r'$\times $'
        ######for y in dummy_y_values:
        ######    yticklabels[np.argwhere(y_array == y)[0][0]] = r'$\times$'


        ###cmap =  sns.color_palette('crest', as_cmap=True).copy()
        cmap =  sns.color_palette('viridis_r', as_cmap=True).copy()
        cmap.set_under('red')
        ####cmap.set_over('black')

        ax = sns.heatmap(np.transpose(z_matrix), linewidth=0.01, xticklabels = yticklabels, yticklabels = xticklabels, cmap = cmap, vmin = vmin, vmax = vmax, norm = LogNorm(vmin = vmin, vmax = vmax))
        ax.set_facecolor('red')

        ## finde index of dummy values to insert separating lines
        x_lines = []
        ix = 0
        for x in dummy_x_values:
            x_lines.append(np.where(x_array_with_dummy == x)[0][0] - ix)
            ix += 1
        y_lines = []
        iy = 0
        for y in dummy_y_values:
            y_lines.append(np.where(y_array_with_dummy == y)[0][0] - iy)
            iy += 1

        ax.hlines(x_lines, *ax.get_xlim(), linewidth = 5, color = 'black')
        ax.vlines(y_lines, *ax.get_ylim(), linewidth = 5, color = 'black')

        ax.invert_yaxis()
        ax.set_xlabel(getPrettyName(var_x))
        ax.set_ylabel(getPrettyName(var_y))


        if "rnorm" in err_type:
            ax.collections[0].colorbar.set_label("Relative error " + getPrettyName(err_type + "_" + self.var_error) )
        else:
            ax.collections[0].colorbar.set_label("Relative error " + getPrettyName(err_type) )

        plt.tight_layout()

        figname = dirname + "/map_" + plot_type + "_" + self.ref_type + "_" + err_type + "_niter" + str(niter);
        if not filter_vars is None:
            figname = setFilenameFilterVars(figname, filter_vars, filter_vals);

        plt.savefig(figname + ".pdf");
        plt.close();

