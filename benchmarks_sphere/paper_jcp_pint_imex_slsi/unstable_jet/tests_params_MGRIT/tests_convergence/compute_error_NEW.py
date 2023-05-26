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
                       "runtime.xbraid_spatial_coarsening",
                       "runtime.xbraid_timestepping_method",
                       "runtime.xbraid_viscosity_coefficient"
                    ]);

####################
## COMPUTE ERRORS ##
####################
es_fine = errorSWEET.errorSWEET(jobs_data, jobs_solution, key_fine, "fine", "sphere", var_error = "prog_phi");
es_fine.computeErrorsPinT(read_errors_from_file = True, pint_type = "xbraid", spectral = True)
es_fine.readPinTResiduals(pint_type = "xbraid");

es_ref = errorSWEET.errorSWEET(jobs_data, jobs_solution, key_fine, "ref", "sphere", var_error = "prog_phi");
es_ref.computeErrorsPinT(read_errors_from_file = True, pint_type = "xbraid", spectral = True)


legend_vars = [
                "runtime.xbraid_max_levels",
                "runtime.xbraid_cfactor",
                ###"runtime.xbraid_nrelax",
                "runtime.xbraid_viscosity_coefficient_coarse",
                "runtime.xbraid_viscosity_coefficient_coarse2",
              ];

groups_vars = [
                "runtime.xbraid_max_levels",
                "runtime.xbraid_cfactor",
                ###"runtime.xbraid_viscosity_coefficient"
                ##"runtime.xbraid_nrelax"
              ];

common_plot_attributes = ["linestyle", "marker"];

filter_vars = [
                "runtime.xbraid_spatial_coarsening",
                "runtime.xbraid_timestepping_method_coarse"
             ];
all_filter_vals = [
                   [51, "l_irk_n_erk"],
                   [128, "l_irk_n_erk"],
                   [51, "l_irk_na_sl_nr_settls_uv_only"],
                   [128, "l_irk_na_sl_nr_settls_uv_only"]
                 ];

remove_vars = [
                 "runtime.xbraid_nrelax",
              ]
remove_vals = [
                 [1, 2, 5],
              ]
es_fine.removeJobs(remove_vars, remove_vals)
es_ref.removeJobs(remove_vars, remove_vals)


## plot errors along iterations for each set of filter vars and for each rnorm
for filter_vals in all_filter_vals:

    for rnorm in [32, 128]:

        ###plot_legend = (rnorm == 32 and 51 in filter_vals)
        ###plot_legend = True

        if "l_irk_n_erk" in filter_vals:
            plot_legend = (rnorm == 32 and 51 in filter_vals)
            ylim = (1e-5, 4e-3)
            legend_fontsize = 16
        else:
            plot_legend = (rnorm == 32 and 128 in filter_vals)
            ylim = (1e-3, 5e-1)
            legend_fontsize = 14

        es_fine.plotErrorPinTAlongIterations(  dirname = "err_Xbraid_fine",
                                               plot_type = "error",
                                               err_type = "rnorm_" + str(rnorm),
                                               legend_vars = legend_vars,
                                               groups_vars = groups_vars,
                                               common_plot_attributes = common_plot_attributes,
                                               filter_vars = filter_vars,
                                               filter_vals = filter_vals,
                                               relative_error = True,
                                               max_iter = 10,
                                               plot_legend = plot_legend,
                                               ncol_legend = 2,
                                               ylim = ylim,
                                               legend_fontsize = legend_fontsize
                                            );

