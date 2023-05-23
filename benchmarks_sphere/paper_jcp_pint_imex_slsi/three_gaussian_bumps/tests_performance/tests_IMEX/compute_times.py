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
                       "runtime.xbraid_timestepping_method",
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

## plot times and speedups in function of error
filter_vars = [
                "runtime.xbraid_timestepping_method_coarse",
                "runtime.xbraid_pt"
             ];
all_filter_vals = [
                   ["l_irk_n_erk", 32]
                 ];

exclude_vars = [
                "runtime.xbraid_timestepping_method_coarse",
                "runtime.xbraid_max_levels",
                "runtime.xbraid_cfactor",
                "runtime.xbraid_nrelax",
                "runtime.xbraid_spatial_coarsening",
                "runtime.xbraid_viscosity_coefficient_coarse1",
                "runtime.xbraid_viscosity_coefficient_coarse2"
             ];
exclude_vals = [
                   ####["l_irk_n_erk", 2, 2, 0, 51, 1e6, "--"],
                   ["l_irk_n_erk", 2, 2, 1, 51, 1e6, "--"],
                   ["l_irk_n_erk", 2, 2, 5, 51, 1e6, "--"],
                   ####["l_irk_n_erk", 2, 2, 0, 128, 1e6, "--"],
                   ["l_irk_n_erk", 2, 2, 1, 128, 1e6, "--"],
                   ####["l_irk_n_erk", 2, 2, 5, 128, 1e6, "--"],
                   ["l_irk_n_erk", 3, 2, 0, 128, 1e27, 1e27],
                   ["l_irk_n_erk", 3, 2, 0, 128, 0, 1e7],
                   ["l_irk_n_erk", 3, 2, 0, 128, 0, 1e17],
                   ["l_irk_n_erk", 3, 2, 0, 128, 0, 1e27],
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

        ts.plotTimesSpeedupsError(   dirname = "times_speedups",
                                     plot_type = plot_type,
                                     err_types = ["rnorm_32", "rnorm_128"],
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
                                     loc_legend1 = 2,
                                     loc_legend2 = 4
                                     ##ylim = (1e-1, 2e1)
                                     ###ylim = (1e-2, 2e1)
                                     ###ncol_legend = 2
                                  );

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
                   ["l_irk_n_erk", 2, 2, 0, 51, 1e6, "--"],
                   ###["l_irk_n_erk", 2, 2, 1, 51, 1e6, "--"],
                   ["l_irk_n_erk", 2, 2, 5, 51, 1e6, "--"],
                   ["l_irk_n_erk", 2, 2, 0, 128, 1e6, "--"],
                   ###["l_irk_n_erk", 2, 2, 1, 128, 1e6, "--"],
                   ["l_irk_n_erk", 2, 2, 5, 128, 1e6, "--"],
                   ###["l_irk_n_erk", 3, 2, 0, 128, 1e27, 1e27],
                   ["l_irk_n_erk", 3, 2, 0, 128, 0, 1e7],
                   ###["l_irk_n_erk", 3, 2, 0, 128, 0, 1e17],
                   ["l_irk_n_erk", 3, 2, 0, 128, 0, 1e27],
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


