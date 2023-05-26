import numpy as np
import matplotlib.pyplot
import sys

###sys.path.append("/home/jcaldass/Development/IME/PostDocIMEUSP/SWEET_scripts/script_paper_JCP");
import dataSWEET
import errorSWEET


t = 60 * 60 * 24 * 6

ds = dataSWEET.dataSWEET("csv", "sphere");
jobs_data, jobs_solution = ds.readAllData("output_prog_phi_t", t_input = t);

ds.storeJobsDetails([
                       "runtime.timestep_size",
                       "runtime.space_res_spectral",
                       "runtime.viscosity",
                       "runtime.viscosity_order"
                    ]);

## find ref simulation (dt = 2)
ref_key = "";
for job in jobs_data.keys():
    if jobs_data[job]["runtime.space_res_spectral"] == 512:
        ref_key = job;
        break;
print("REF:", ref_key);


es = errorSWEET.errorSWEET(jobs_data, jobs_solution, ref_key, "ref", "sphere", var_error = "prog_phi");
es.computePhysicalErrors(t = t)
es.computeSpectralErrors(t = t)

filter_vars = [
              ]

filter_vals = [
              ]

legend_vars = [
                "runtime.space_res_spectral",
                "runtime.timestep_size"
              ]

legend_vars_plot_attributes = [
                                 ["marker"],                               ## same res_spectral -> same marker and linestyle
                                 ["color", "linestyle"]                    ## same color -> same color
                              ]

for err_type in ["err_L1", "err_L2", "err_Linf", "rnorm_16", "rnorm_32", "rnorm_64", "rnorm_128"]:
    es.plotErrorSerial(  dirname = "err_serial",
                         err_type = err_type,
                         x_var = "runtime.timestep_size",
                         legend_vars = None,
                         groups_vars = None,
                         filter_vars = filter_vars,
                         filter_vals = filter_vals,
                         legend_vars_plot_attributes = legend_vars_plot_attributes
                      )

es.plotErrorSerial(  dirname = "err_serial",
                     err_type = ["err_L2", "rnorm_16", "rnorm_32", "rnorm_64", "rnorm_128"],
                     x_var = "runtime.timestep_size",
                     legend_vars = None,
                     groups_vars = None,
                     filter_vars = filter_vars,
                     filter_vals = filter_vals,
                     legend_vars_plot_attributes = legend_vars_plot_attributes
                  )


