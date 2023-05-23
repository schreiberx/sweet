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
import spectrumSWEET

###############
## READ DATA ##
###############
ds = dataSWEET.dataSWEET("csv", "sphere");

jobs_data, jobs_solution = ds.readAllData("output_prog_phi_pert_t");

key_fine, key_ref = ds.findFineSimulation("xbraid");
print("\tfine simulation: ", key_fine, jobs_data[key_fine]["jobgeneration.p_job_dirpath"]);
print("\tref simulation: ", key_ref, jobs_data[key_fine]["jobgeneration.p_job_dirpath"]);

print("Number of jobs: ", len(ds.getListJobs()));

ds.storeJobsDetails([
                       "runtime.timestep_size",
                       "runtime.xbraid_enabled",
                       "runtime.xbraid_max_levels",
                       "runtime.xbraid_cfactor",
                       "runtime.xbraid_nrelax",
                       "runtime.xbraid_spatial_coarsening",
                       "runtime.xbraid_timestepping_method"
                    ]);

###############################
## READ AND PLOT KE SPECTRUM ##
###############################
ss = spectrumSWEET.spectrumSWEET(jobs_data, "sphere", key_fine, key_ref)

legend_vars = [
                "runtime.xbraid_max_levels",
                "runtime.xbraid_cfactor",
                "runtime.xbraid_nrelax",
                "runtime.xbraid_spatial_coarsening"
              ];

niters = [0, 5, 7]
niters = [0, 5, 7]
ss.computeKESpectrum(t = 60. * 60. * 36., niters = niters)
ss.plotKESpectrum(title = "", dirname = "spectrum_plots", output_filename = "spec", legend_vars = legend_vars, PinT = True, niters = niters)

