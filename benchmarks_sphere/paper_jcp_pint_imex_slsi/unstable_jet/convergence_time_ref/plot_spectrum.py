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

###sys.path.append("/home/jcaldass/Development/IME/PostDocIMEUSP/SWEET_scripts/script_paper_JCP")
import dataSWEET
import spectrumSWEET

###############
## READ DATA ##
###############
ds = dataSWEET.dataSWEET("csv", "sphere");

jobs_data, jobs_solution = ds.readAllData("output_prog_phi_pert_t");


## find ref simulation (dt = 2)
ref_key = "";
for job in jobs_data.keys():
    if "job_benchref" in jobs_data[job]["runtime.p_job_dirpath"]:
        ref_key = job;
        break;
print("REF:", ref_key);

print("Number of jobs: ", len(ds.getListJobs()));

ds.storeJobsDetails([
                       "runtime.timestep_size",
                       "runtime.space_res_spectral"
                    ]);

###############################
## READ AND PLOT KE SPECTRUM ##
###############################
ss = spectrumSWEET.spectrumSWEET(jobs_data, "sphere", ref_key, ref_key)

legend_vars = [
                "runtime.timestep_size",
                "runtime.space_res_spectral"
              ];

ss.computeKESpectrum(t = 24. * 6. * 60 * 60, niters = None)
ss.plotKESpectrum(title = "", dirname = "spectrum_plots", output_filename = "spec", legend_vars = legend_vars, PinT = False, niters = [0])

