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


class timesSWEET:

    def __init__(self, jobs_data, jobs_errors, ref_key, ref_type, geometry, var_error = "prog_phi_pert"):

        self.jobs_data = jobs_data
        self.jobs_errors = jobs_errors
        self.ref_key = ref_key
        self.ref_type = ref_type
        self.geometry = geometry
        self.var_error = var_error

        self.times = {};
        self.speedups = {};
        self.times_error = {};
        self.speedups_error = {};

    ## some jobs were run to compute times
    ## and others to compute errors
    ## this function makes the link between jobs using same parameters
    def linkErrorAndTimeJobs(self, list_vars):

        ## loop over time jobs
        for job1 in self.jobs_data.keys():

            if job1 == self.ref_key:
                continue

            ## jump jobs containing output (errors)
            if not ("runtime.xbraid_no_output" in self.jobs_data[job1]):
                continue
            if not self.jobs_data[job1]["runtime.xbraid_no_output"]:
                continue

            ## find job with same parameters
            for job2 in self.jobs_data.keys():

                if job2 == self.ref_key:
                    continue

                ## skip time jobs
                if "runtime.xbraid_no_output" in self.jobs_data[job2]:
                    if self.jobs_data[job2]["runtime.xbraid_no_output"]:
                        continue

                ###print("\n\nComparing " + getJobPath(self.jobs_data, job1) + " " + getJobPath(self.jobs_data, job2))
                job_found = True
                for var in list_vars:
                    ###print(var, self.jobs_data[job1][var], self.jobs_data[job2][var])
                    if not (self.jobs_data[job1][var] == self.jobs_data[job2][var]):
                        job_found = False
                        break

                if job_found:
                    break

            if not job_found:
                raise Exception("Not found corresponding job to " + getJobPath(self.jobs_data, job1))

            ## copy errors to time job
            assert not (job1 in self.jobs_errors.err_iters_array.keys())
            self.jobs_errors.err_iters_array[job1] = self.jobs_errors.err_iters_array[job2].copy()

    ## eclude jobs run to compute errors (and not computing times)
    def excludeErrorJobs(self):

        print("Nb of jobs before exclude:", len(self.jobs_data) )
        for job in self.jobs_data.copy().keys():
            if job == self.ref_key:
                continue
            exclude_job = False
            if not ("runtime.xbraid_no_output" in self.jobs_data[job].keys()):
                exclude_job = True
            elif not self.jobs_data[job]["runtime.xbraid_no_output"]:
                exclude_job = True

            if exclude_job:
                if job in self.jobs_data.keys():
                    del self.jobs_data[job]
        print("Nb of jobs after exclude:", len(self.jobs_data) )

    ## create file with jobs details (values of given list of parameters)
    def storeJobsDetails(self, list_params):

        f = open("job_details_time", "w");
        print("Storing details of " + str(len(self.jobs_data)) + " jobs")
        for job in self.jobs_data.keys():
            f.write("\n" + getJobPath(self.jobs_data, job) + "\n");
            for param in list_params:
                try:
                    val = self.jobs_data[job][param];
                    f.write(param + "\t\t\t\t\t" + str(val) + "\n");
                except:
                    pass;
            f.write("\n");
        f.close();


    def readPinTTimes(self):

        for job in self.jobs_data.keys():

            if job == self.ref_key:
                continue

            fname = getJobPath(self.jobs_data, job) + "/output.out"

            if not os.path.isfile(fname):
                print(getJobPath(self.jobs_data, job))
                self.times[job] = None
                continue

            self.times[job] = np.empty((0,2), dtype = float);

            ## read errors
            lines = [line.rstrip() for line in open(fname)];

            for line in lines:
                spl = line.split()
                if "Braid:" in spl and "wall" in spl and "time" in spl:
                    niter = int(spl[2][2:]);
                    time = float(spl[-1])
                    if niter in self.times[job][:, 0]:
                        self.times[job][niter][1] = time
                    else:
                        self.times[job] = np.vstack((self.times[job], np.array([niter, time])))

    def readRefTime(self):

        fname = getJobPath(self.jobs_data, self.ref_key) + "/output.out"

        lines = [line.rstrip() for line in open(fname)];

        for line in lines:
            spl = line.split()

            if "simulation_benchmark_timings.main:" in spl:
                self.times[self.ref_key] = float(spl[-1]);
                break;

    def computeSpeedup(self):

        for job in self.times.keys():

            if job == self.ref_key:
                continue

            self.speedups[job] = np.copy(self.times[job])
            self.speedups[job][:, 1] =  self.times[self.ref_key] / self.speedups[job][:, 1]

    def computeTimesSpeedupsInFunctionErrors(self):

        base_err_thresholds = np.arange(0, -9, -1, dtype = float);
        base_err_thresholds = np.array([10 ** t for t in base_err_thresholds])
        err_thresholds = np.append(base_err_thresholds, np.round(.5 * base_err_thresholds, 10))
        err_thresholds = np.append(err_thresholds, np.round(.25 * base_err_thresholds, 10))
        err_thresholds = np.append(err_thresholds, np.round(.75 * base_err_thresholds, 10))
        ###err_thresholds = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
        ###err_types = ["err_L1", "err_L2", "err_Linf", "rnorm_16", "rnorm_32", "rnorm_64", "rnorm_128", "rnorm_256"]

        for job in self.times.keys():

            if job == self.ref_key:
                continue

            self.times_error[job] = {}
            self.speedups_error[job] = {}

            all_errors = self.jobs_errors.err_iters_array[job]
            for err_type in all_errors.keys():

                ## get errors
                errors = all_errors[err_type]

                ###self.times_error[job][err_type] = np.empty((0, 2), dtype = float)
                ###self.speedups_error[job][err_type] = np.empty((0, 2), dtype = float)
                self.times_error[job][err_type] = {}
                self.speedups_error[job][err_type] = {}

                for err_threshold in err_thresholds:

                    niter_threshold = 1e16;
                    for niter in range(errors.shape[0]):
                        if errors[niter, 1] < err_threshold and niter < niter_threshold:
                            niter_threshold = niter;
                    ###print(self.times_error[job][err_type])
                    ###print(niter_threshold)
                    if niter_threshold < errors.shape[0]:
                        time = self.times[job][niter_threshold, 1]
                        speedup = self.speedups[job][niter_threshold, 1]
                    else:
                        time = 1e16
                        speedup = 1e-16

                    ###self.times_error[job][err_type] = np.vstack((self.times_error[job][err_type], np.array([err_threshold, time])))
                    ###self.speedups_error[job][err_type] = np.vstack((self.speedups_error[job][err_type], np.array([err_threshold, speedup])))
                    self.times_error[job][err_type][err_threshold] = time
                    self.speedups_error[job][err_type][err_threshold] = speedup

                ####print(self.speedups_error[job][err_type])

    def getTimesSpeedup(self):

        self.readRefTime();
        self.readPinTTimes();
        self.computeSpeedup();
        ####self.computeTimesSpeedupsInFunctionErrors()

        print("Nb of jobs with computed time: ", len(self.times), len(self.speedups))

    def plotTimesSpeedupsAlongIterations(self, dirname, plot_type, legend_vars, groups_vars, common_plot_attributes, filter_vars = None, filter_vals = None, max_iter = None, plot_legend = True, ncol_legend = 1, first_plot_idx = 0, ylim = None, linewidth = .75):

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        if plot_type == "times":
            array_jobs = self.times.copy();
        elif plot_type == "speedups":
            array_jobs = self.speedups.copy();
        else:
            sys.exit("Wrong plot_type: " + plot_type);

        if self.ref_key in array_jobs:
            del array_jobs[self.ref_key]

        jobs_to_plot = sortJobs(array_jobs, legend_vars, self.jobs_data)

        jobs_to_plot = filterJobs(jobs_to_plot, filter_vars, filter_vals, self.jobs_data)

        ## group jobs based on given parameters
        groups = createGroups(jobs_to_plot, groups_vars, self.jobs_data);

        fig, ax = plt.subplots();
        igroup = first_plot_idx;
        for group in groups.values():

            ijob = 0;
            for job in group:

                array = array_jobs[job];

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

                label = getLabel(legend_vars, job, self.jobs_data)
                ax.plot(array[:range_iter, 0], array[:range_iter, 1], label = label, color = cml["color"], marker = cml["marker"], linestyle = cml["linestyle"], linewidth = linewidth)
                ijob += 1

            igroup += 1;

        horizontal_line_x = np.arange(0, max_iter + 1, 1)
        if plot_type == "times":
            horizontal_line_y = np.ones_like(horizontal_line_x) * self.times[self.ref_key]
        else:
            horizontal_line_y = np.ones_like(horizontal_line_x)

        ax.plot(horizontal_line_x, horizontal_line_y, linestyle = "--", color = 'red', label = "ref")

        if not (ylim is None):
            ax.set_ylim(ylim)

        ax.set_yscale("log");
        ax.set_xlabel("Iteration");
        if plot_type == "times":
            ax.set_ylabel("Wall times (s)" );
        else:
            ax.set_ylabel("Speedup");
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        if plot_legend:
            ax.legend(title = setTitleLegend(legend_vars), ncol = ncol_legend, loc='best');

        plt.tight_layout();

        figname = dirname + "/" + plot_type + "_" + self.ref_type;
        if not filter_vars is None:
            figname = setFilenameFilterVars(figname, filter_vars, filter_vals);

        plt.savefig(figname + ".pdf");
        plt.close();

    def plotTimesSpeedupsInFunctionProcessors(self, dirname, plot_type, common_plot_attributes, niters, filter_vars = None, filter_vals = None, ylim = None, linewidth = .75):

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        if plot_type == "times":
            array_jobs = self.times.copy();
        elif plot_type == "speedups":
            array_jobs = self.speedups.copy();
        else:
            sys.exit("Wrong plot_type: " + plot_type);

        if self.ref_key in array_jobs:
            del array_jobs[self.ref_key]

        jobs_to_plot = filterJobs(list(array_jobs.keys()), filter_vars, filter_vals, self.jobs_data)


        max_nb_pt = 0;

        fig, ax = plt.subplots();

        ## collect numbers of processors and errors
        initer = 0
        for niter in niters:

            ijob = 0;

            plot_x = np.empty(0, dtype = int)
            plot_y = np.empty(0, dtype = float)

            for job in jobs_to_plot:

                array = array_jobs[job];

                nb_pt = self.jobs_data[job]["runtime.xbraid_pt"];
                print(getJobPath(self.jobs_data, job))
                y = array_jobs[job][niter, 1]

                plot_x = np.append(plot_x, nb_pt)
                plot_y = np.append(plot_y, y)

                ijob += 1

            ## sort
            idx = np.argsort(plot_x)
            plot_x = plot_x[idx]
            plot_y = plot_y[idx]

            ## get color, marker, linestyle
            cml = {};
            for pa in plot_attributes.keys():
                cml[pa] = getPlotAttribute(pa, initer);

            ## plot
            label = r'$k = {}$'.format(niter)
            ax.plot(plot_x, plot_y, label = label, color = cml["color"], marker = cml["marker"], linestyle = cml["linestyle"], linewidth = linewidth)

            max_nb_pt = max(max_nb_pt, np.max(plot_x))

            initer += 1

        if (max_nb_pt < 100):
            ax.set_xlim(.75, 100)
            max_nb_pt = 99

        horizontal_line_x = np.arange(0, max_nb_pt + 1, 1)
        if plot_type == "times":
            horizontal_line_y = np.ones_like(horizontal_line_x) * self.times[self.ref_key]
        else:
            horizontal_line_y = np.ones_like(horizontal_line_x)

        ax.plot(horizontal_line_x, horizontal_line_y, linestyle = "--", color = 'red')

        if not (ylim is None):
            ax.set_ylim(ylim)


        ax.set_xscale("log");
        ax.set_yscale("log");
        ax.set_xlabel(getPrettyName('runtime.xbraid_pt'));
        if plot_type == "times":
            ax.set_ylabel("Wall times (s)" );
        else:
            ax.set_ylabel("Speedup");
        ###ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        ax.legend();

        plt.tight_layout();

        figname = dirname + "/" + plot_type + "_" + self.ref_type + "_nbpt";
        if not filter_vars is None:
            figname = setFilenameFilterVars(figname, filter_vars, filter_vals);

        plt.savefig(figname + ".pdf");
        plt.close();



    def plotTimesSpeedupsError(self, dirname, plot_type, err_types, legend_vars, groups_vars, common_plot_attributes, err_thresholds = None, niters = None, filter_vars = None, filter_vals = None, exclude_vars = None, exclude_vals = None, max_iter = None, plot_legend = True, ncol_legend = 1, first_plot_idx = 0, xlim = None, ylim = None, linewidth = .75, loc_legend1 = 'best', loc_legend2 = 'best'):

        assert not(err_thresholds is None and niters is None)

        if not os.path.isdir(dirname):
            os.makedirs(dirname);

        if not (err_thresholds is None):
            if plot_type == "times":
                array_jobs = self.times_error.copy();
            elif plot_type == "speedups":
                array_jobs = self.speedups_error.copy();
            else:
                sys.exit("Wrong plot_type: " + plot_type);
        elif not (niters is None):
            if plot_type == "times":
                array_jobs = self.times.copy();
            elif plot_type == "speedups":
                array_jobs = self.speedups.copy();
            else:
                sys.exit("Wrong plot_type: " + plot_type);

        if self.ref_key in array_jobs:
            del array_jobs[self.ref_key]

        jobs_to_plot = sortJobs(array_jobs, legend_vars, self.jobs_data)

        jobs_to_plot = filterJobs(jobs_to_plot, filter_vars, filter_vals, self.jobs_data)

        jobs_to_plot = excludeJobs(jobs_to_plot, exclude_vars, exclude_vals, self.jobs_data)

        ## group jobs based on given parameters
        groups = createGroups(jobs_to_plot, groups_vars, self.jobs_data);

        min_x = 1e10
        max_x = -1e10
        fig, ax = plt.subplots();
        igroup = first_plot_idx;
        for group in groups.values():

            ijob = 0;
            for job in group:

                ierr_type = 0;
                for err_type in err_types:

                    x_plot = np.empty(0, dtype = float)
                    y_plot = np.empty(0, dtype = float)

                    if not (err_thresholds is None):
                        array = array_jobs[job][err_type];
                    elif not (niters is None):
                        array = array_jobs[job];

                    ## get color, marker, linestyle
                    cml = {};
                    for pa in plot_attributes.keys():
                        if pa in common_plot_attributes:
                            cml[pa] = getPlotAttribute(pa, igroup);
                        else:
                            cml[pa] = getPlotAttribute(pa, ierr_type);

                    if not (err_thresholds is None):
                        for err_threshold in err_thresholds:
                            val = array[err_threshold]
                            if (plot_type == "times" and val > 1e10) or (plot_type == "speedups" and val < 1e-10):
                                continue
                            x_plot = np.append(x_plot, err_threshold)
                            y_plot = np.append(y_plot, array[err_threshold])

                    elif not (niters is None):
                         if array.shape[0] == 0:
                             continue
                         for niter in niters:
                            if niter >= array.shape[0]:
                                break
                            if niter >= self.jobs_errors.err_iters_array[job][err_type].shape[0]:
                                break
                            ###print(self.jobs_errors.err_iters_array[job][err_type], getJobPath(self.jobs_data, job))
                            ###print(array, getJobPath(self.jobs_data, job))
                            x_plot = np.append(x_plot, self.jobs_errors.err_iters_array[job][err_type][niter, 1])
                            y_plot = np.append(y_plot, array[niter, 1])

                         min_x = np.min([min_x, min(x_plot)])
                         max_x = np.max([max_x, max(x_plot)])

                    if ierr_type == 0:
                        label = getLabel(legend_vars, job, self.jobs_data)
                    else:
                        label = ""
                    ax.plot(x_plot, y_plot, label = label, color = cml["color"], marker = cml["marker"], linestyle = cml["linestyle"], linewidth = linewidth)

                    ierr_type += 1

                ijob += 1

            igroup += 1;


        if not (xlim is None):
            horizontal_line_x = np.array([xlim[0], xlim[1]])
        else:
            if not (err_thresholds is None):
                horizontal_line_x = np.array([min(err_thresholds), max(err_thresholds)])
            elif not (niters is None):
                horizontal_line_x = np.array([min_x, max_x * 1.5])

        if plot_type == "times":
            horizontal_line_y = np.ones_like(horizontal_line_x) * self.times[self.ref_key]
        else:
            horizontal_line_y = np.ones_like(horizontal_line_x)

        ax.plot(horizontal_line_x, horizontal_line_y, linestyle = "-.", color = 'red', linewidth = linewidth)

        if not (xlim is None):
            ax.set_xlim(xlim)
        if not (ylim is None):
            ax.set_ylim(ylim)

        ax.set_xscale("log");
        ax.set_yscale("log");
        err_threshold_str = ""
        if not (err_thresholds is None):
            err_threshold_str = " threshold"
        if len(err_types) == 1:
            if "rnorm" in err_types[0]:
                ax.set_xlabel("Relative error " + getPrettyName(err_types[0] + "_" + self.var_error) + err_threshold_str);
            else:
                ax.set_xlabel("Relative error " + getPrettyName(err_types[0]) + err_threshold_str);
        else:
            if "rnorm" in err_types[0]:
                ax.set_xlabel("Relative error " + getPrettyName("rnorm_" + self.var_error) + err_threshold_str);
            else:
                ax.set_xlabel("Relative error" + err_threshold_str);

        if plot_type == "times":
            ax.set_ylabel("Wall times (s)" );
        else:
            ax.set_ylabel("Speedup");

        if plot_legend:
            legend = ax.legend(title = setTitleLegend(legend_vars), ncol = ncol_legend, loc = loc_legend1, fontsize  = 14, title_fontsize = 12);
            ###ax.legend(ncol = ncol_legend, loc='best', fontsize  = 14);
            ###ax.legend(title = setTitleLegend(legend_vars), ncol = ncol_legend, loc='upper center', bbox_to_anchor=(0.5, 1.05) );

            ## add second legend if there is more than one err_type
            if len(err_types) > 1:
                styles2 = []
                labels2 = []

                ierr_type = 0
                for err_type in err_types:

                    cml = {};
                    cml["marker"] = "None"
                    cml["color"] = "black"
                    cml["linestyle"] = "--"
                    for pa in plot_attributes.keys():
                        if not (pa in common_plot_attributes):
                            cml[pa] = getPlotAttribute(pa, ierr_type);

                    styles2.append(plt.Line2D((0, 1), (0, 0), color = cml["color"], marker = cml["marker"], linestyle = cml["linestyle"]));
                    labels2.append(r'$R_{\mathrm{norm}} = $' + " " + r'${}$'.format(int(err_type[6:])))

                    ierr_type += 1

                legend2 = plt.legend(styles2, labels2, loc = loc_legend2, fontsize = 14);
                ax.add_artist(legend)
                ax.add_artist(legend2);




        plt.tight_layout();

        if not (err_thresholds is None):
            figname = dirname + "/errthresholds_" + plot_type + "_" + self.ref_type + "_" + err_type[0];
        elif not (niters is None):
            figname = dirname + "/errs_" + plot_type + "_" + self.ref_type + "_" + err_types[0];
        if not filter_vars is None:
            figname = setFilenameFilterVars(figname, filter_vars, filter_vals);

        plt.savefig(figname + ".pdf");
        plt.close();

