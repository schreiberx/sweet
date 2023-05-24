import numpy as np
import sys
import math
import scipy
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
from matplotlib.transforms import Transform
from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator)
from matplotlib.lines import Line2D
import os

####sys.path.append("/home/jcaldass/Development/IME/PostDocIMEUSP/SWEET_scripts/script_paper_JCP")
import dataSWEET
import errorSWEET
import timesSWEET

###############
## READ DATA ##
###############
ds = dataSWEET.dataSWEET("csv", "sphere");
####ds.read_data = False;
####ds.read_xbraid_parameters = True;

jobs_data, jobs_solution = ds.readAllData("output_prog_phi_t");

jobs_data = ds.splitPinTParameters("xbraid")

key_fine, key_ref = ds.findFineSimulation("xbraid");
print("\tfine_simulation: ", key_fine, jobs_data[key_fine]["jobgeneration.p_job_dirpath"]);

print("Number of jobs: ", len(ds.getListJobs()));

ds.storeJobsDetails([
                       "runtime.timestep_size",
                       "runtime.xbraid_enabled",
                       "runtime.xbraid_max_levels",
                       "runtime.xbraid_cfactor",
                       "runtime.xbraid_nrelax",
                       "runtime.xbraid_spatial_coarsening"
                    ]);

####################
## COMPUTE ERRORS ##
####################
es_fine = errorSWEET.errorSWEET(jobs_data, jobs_solution, key_fine, "fine", "sphere", var_error = "prog_phi");
es_fine.computeErrorsPinT(read_errors_from_file = True, pint_type = "xbraid", spectral = True)
es_fine.readPinTResiduals(pint_type = "xbraid");

es_ref = errorSWEET.errorSWEET(jobs_data, jobs_solution, key_fine, "ref", "sphere", var_error = "prog_phi");
es_ref.computeErrorsPinT(read_errors_from_file = True, pint_type = "xbraid", spectral = True)


###############
## GET TIMES ##
###############
ts = timesSWEET.timesSWEET(jobs_data, es_fine, key_fine, "fine", "sphere", var_error = "prog_phi");

vars_common_jobs = {
                       "runtime.timestep_size",
                       "runtime.xbraid_enabled",
                       "runtime.xbraid_max_levels",
                       "runtime.xbraid_cfactor",
                       "runtime.xbraid_nrelax",
                       "runtime.xbraid_spatial_coarsening",
                       "runtime.xbraid_viscosity_order_fine",
                       "runtime.xbraid_viscosity_order_coarse",
                       "runtime.xbraid_viscosity_order_coarse2",
                       "runtime.xbraid_viscosity_coefficient_fine",
                       "runtime.xbraid_viscosity_coefficient_coarse",
                       "runtime.xbraid_viscosity_coefficient_coarse2"
                   }
ts.linkErrorAndTimeJobs(list_vars = vars_common_jobs)
ts.excludeErrorJobs()
ts.getTimesSpeedup()
ts.storeJobsDetails([
                       "runtime.xbraid_max_levels",
                       "runtime.xbraid_cfactor",
                       "runtime.xbraid_nrelax",
                       "runtime.xbraid_spatial_coarsening",
                        "runtime.xbraid_viscosity_coefficient_coarse1",
                        "runtime.xbraid_viscosity_coefficient_coarse2",
                        "runtime.xbraid_pt",
                    ]);


common_plot_attributes = ["color", "marker"];

## plot times and speedups along iterations
filter_vars = [
                "runtime.xbraid_timestepping_method_coarse",
                "runtime.xbraid_pt"
             ];
all_filter_vals = [
                   ["l_irk_na_sl_nr_settls_uv_only", 1],
                   ["l_irk_na_sl_nr_settls_uv_only", 2],
                   ["l_irk_na_sl_nr_settls_uv_only", 4],
                   ["l_irk_na_sl_nr_settls_uv_only", 8],
                   ["l_irk_na_sl_nr_settls_uv_only", 16],
                 ];

for filter_vals in all_filter_vals:

    for plot_type in ["times", "speedups"]:

        plot_legend = True

        legend_vars = [
                        "runtime.xbraid_max_levels",
                        "runtime.xbraid_cfactor",
                        "runtime.xbraid_nrelax",
                        "runtime.xbraid_spatial_coarsening",
                        "runtime.xbraid_viscosity_coefficient_coarse1",
                        "runtime.xbraid_viscosity_coefficient_coarse2",
                      ];

        ts.plotTimesSpeedupsAlongIterations(   dirname = "times_speedups",
                                               plot_type = plot_type,
                                               legend_vars = legend_vars,
                                               groups_vars = legend_vars,
                                               common_plot_attributes = common_plot_attributes,
                                               filter_vars = filter_vars,
                                               filter_vals = filter_vals,
                                               max_iter = 10,
                                               plot_legend = plot_legend,
                                               ncol_legend = 2
                                               ###ncol_legend = 2
                                            );


## plot times and speedups in function of error thresholds
filter_vars = [
                "runtime.xbraid_timestepping_method_coarse",
                "runtime.xbraid_pt"
             ];
all_filter_vals = [
                   ##["l_irk_na_sl_nr_settls_uv_only", 16]
                   ["l_irk_na_sl_nr_settls_uv_only", 8]
                 ];

exclude_vars = [
             ];
exclude_vals = [
                 ];


for filter_vals in all_filter_vals:

    for plot_type in ["times", "speedups"]:

        for err_type in ["rnorm_32", "rnorm_128"]:

            plot_legend = (err_type == "rnorm_32")

            legend_vars = [
                            "runtime.xbraid_max_levels",
                            "runtime.xbraid_cfactor",
                            "runtime.xbraid_nrelax",
                            "runtime.xbraid_spatial_coarsening",
                            "runtime.xbraid_viscosity_coefficient_coarse1",
                            "runtime.xbraid_viscosity_coefficient_coarse2",
                          ];

            ts.plotTimesSpeedupsErrorThresholds(   dirname = "times_speedups",
                                                   plot_type = plot_type,
                                                   err_type = err_type,
                                                   ##err_thresholds = [1e-2, 7.5e-3, 5e-3, 2.5e-3, 1e-3, 7.5e-4, 5e-4, 2.5e-4, 1e-4, 7.5e-5, 5e-5, 2.5e-5],
                                                   niters = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                                   legend_vars = legend_vars,
                                                   groups_vars = legend_vars,
                                                   common_plot_attributes = common_plot_attributes,
                                                   filter_vars = filter_vars,
                                                   filter_vals = filter_vals,
                                                   exclude_vars = exclude_vars,
                                                   exclude_vals = exclude_vals,
                                                   max_iter = 10,
                                                   plot_legend = plot_legend,
                                                   ncol_legend = 1,
                                                   ###ylim = (1e-1, 2e1)
                                                   ###ylim = (1e-2, 2e1)
                                                   ###ncol_legend = 2
                                                );


    ###def plotTimesSpeedupsErrorThreshold(self, dirname, plot_type, err_type, err_thresholds, legend_vars, groups_vars, common_plot_attributes, filter_vars = None, filter_vals = None, max_iter = None, plot_legend = True, ncol_legend = 1, first_plot_idx = 0, ylim = None, linewidth = .75):


## plot times and speedups in function of processors
filter_vars = [
                "runtime.xbraid_timestepping_method_coarse",
                "runtime.xbraid_max_levels",
                "runtime.xbraid_cfactor",
                "runtime.xbraid_nrelax",
                "runtime.xbraid_spatial_coarsening",
                "runtime.xbraid_viscosity_coefficient_coarse1",
                "runtime.xbraid_viscosity_coefficient_coarse2",

             ];
all_filter_vals = [
                   ["l_irk_na_sl_nr_settls_uv_only", 2, 2, 0, 51, 1e6, "--"],
                  ];

for filter_vals in all_filter_vals:

    for plot_type in ["times", "speedups"]:

        plot_legend = True

        legend_vars = [
                        "runtime.xbraid_max_levels"
                      ];

        ts.plotTimesSpeedupsInFunctionProcessors(   dirname = "times_speedups",
                                                    plot_type = plot_type,
                                                    common_plot_attributes = common_plot_attributes,
                                                    niters = [0, 1, 5],
                                                    filter_vars = filter_vars,
                                                    filter_vals = filter_vals
                                                 );


######plotTimesSpeedupsInFunctionProcessors(self, dirname, plot_type, legend_vars, groups_vars, common_plot_attributes, niter, filter_vars = None, filter_vals = None, plot_legend = True, ncol_legend = 1, first_plot_idx = 0, ylim = None, linewidth = .75):
####
####legend_vars = [
####                "runtime.xbraid_max_levels",
####                "runtime.xbraid_cfactor",
####                "runtime.xbraid_nrelax"
####              ];
####
####groups_vars = [
####                "runtime.xbraid_max_levels",
####                "runtime.xbraid_cfactor",
####                "runtime.xbraid_nrelax",
####                "runtime.xbraid_spatial_coarsening"
####              ];
####
####common_plot_attributes = ["linestyle", "marker"];
####
####filter_vars = [
####                "runtime.xbraid_timestepping_method_coarse"
####
####             ];
####all_filter_vals = [
####                   ["l_irk_n_erk"]
####                 ];
####
####remove_vars = [
####                 "runtime.xbraid_nrelax"
####              ]
####remove_vals = [
####                 [2]
####              ]
####es_fine.removeJobs(remove_vars, remove_vals)
####es_ref.removeJobs(remove_vars, remove_vals)
####
####
####
####
###### plot errors along iterations for each set of filter vars and for each rnorm
####for filter_vals in all_filter_vals:
####
####    ###for rnorm in [16, 32, 64, 128, 256]:
####    for rnorm in [32, 128]:
####
####        plot_legend = (rnorm == 32 and 51 in filter_vals)
####
####        es_fine.plotErrorPinTAlongIterations(  dirname = "err_Xbraid_fine",
####                                               plot_type = "error",
####                                               err_type = "rnorm_" + str(rnorm),
####                                               legend_vars = legend_vars,
####                                               groups_vars = groups_vars,
####                                               common_plot_attributes = common_plot_attributes,
####                                               filter_vars = filter_vars,
####                                               filter_vals = filter_vals,
####                                               relative_error = True,
####                                               max_iter = 10,
####                                               plot_legend = plot_legend,
####                                               ncol_legend = 2
####                                               ###ncol_legend = 2
####                                            );
####
####        es_ref.plotErrorPinTAlongIterations(   dirname = "err_Xbraid_ref",
####                                               plot_type = "error",
####                                               err_type = "rnorm_" + str(rnorm),
####                                               legend_vars = legend_vars,
####                                               groups_vars = groups_vars,
####                                               common_plot_attributes = common_plot_attributes,
####                                               filter_vars = filter_vars,
####                                               filter_vals = filter_vals,
####                                               relative_error = True,
####                                               max_iter = 10,
####                                               plot_legend = plot_legend,
####                                               ncol_legend = 2
####                                               ###ncol_legend = 2
####                                            );
