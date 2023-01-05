#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob

from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *


def read_error_file(path, error_type):

    lines = [line.rstrip() for line in open(path)];

    if error_type == "physical":

        err_L1 = -1;
        err_L2 = -1;
        err_Linf = -1;

        for line in lines:
            spl = line.split();
            if spl[0] == "errL1":
                err_L1 = float(spl[1]);
                continue;
            if spl[0] == "errL2":
                err_L2 = float(spl[1]);
                continue;
            if spl[0] == "errLinf":
                err_Linf = float(spl[1]);
                continue;

        if err_L1 < 0:
            raise Exception("ERROR: err_L1 not found in " + path);
        if err_L2 < 0:
            raise Exception("ERROR: err_L2 not found in " + path);
        if err_Linf < 0:
            raise Exception("ERROR: err_Linf not found in " + path);

        return err_L1, err_L2, err_Linf

    elif error_type == "spectral":
        err = {};

        for line in lines:
            spl = line.split();
            if spl[0] == "errLinf":
                err[int(spl[1])] = float(spl[2]);

        return err;

    else:
        raise Exception("Wrong error type")


path_simulations = sys.argv[1];
tmp_fine_sim = sys.argv[2];

## get list of jobs in this directory
list_jobs = glob.glob(path_simulations + "/job_bench_*");
list_jobs = [os.path.basename(job) for job in list_jobs];
## exclude fine simulation
if tmp_fine_sim in list_jobs:
    list_jobs.remove(tmp_fine_sim);

print("    ** {} jobs found.".format(len(list_jobs)));


list_vars_parareal = [
                        "runtime.timestep_size", "runtime.timestepping_order", "runtime.timestepping_order2", "runtime.max_simulation_time",
                        "runtime.parareal_enabled", "runtime.parareal_coarse_slices", "runtime.parareal_coarse_timestepping_method", "runtime.parareal_coarse_timestepping_order"
                     ];
list_vars_xbraid = [
                        "runtime.timestep_size", "runtime.timestepping_order", "runtime.timestepping_order2", "runtime.max_simulation_time",
                        "runtime.xbraid_enabled", "runtime.xbraid_cfactor", "runtime.xbraid_max_levels", "runtime.xbraid_pt", "runtime.xbraid_store_iterations",
                     ];
list_vars_common = [
                        "runtime.timestep_size", "runtime.timestepping_order", "runtime.timestepping_order2", "runtime.max_simulation_time"
                   ];


jd = JobsData(path_simulations + '/job_bench_*', verbosity=0).get_flattened_data();
## get useful job info
job_info = {};
for key in jd.keys():
    path = os.path.basename(jd[key]["jobgeneration.p_job_dirpath"]);
    if path == tmp_fine_sim:
        continue;
    job_info[path] = {};
    if jd[key]["runtime.parareal_enabled"]:
        for s in list_vars_parareal:
            job_info[path][s] = jd[key][s];
        job_info[path]["runtime.xbraid_enabled"] = False;
    elif jd[key]["runtime.xbraid_enabled"]:
        for s in list_vars_xbraid:
            job_info[path][s] = jd[key][s];
        job_info[path]["runtime.parareal_enabled"] = False;
    else:
        sys.exit("Simulation " + path + "must be parareal or xbraid.");



## find identical jobs
small = 1e-8;
read_jobs = [];
ipair = 0;
for job1 in list_jobs:
    if job1 in read_jobs:
        continue;
    read_jobs.append(job1);
    if job_info[job1]["runtime.parareal_enabled"]:
        p_or_x = "parareal";
    elif job_info[job1]["runtime.xbraid_enabled"]:
        p_or_x = "xbraid";

    found_job = False;
    for job2 in list_jobs:
        if job2 in read_jobs:
            continue;
        found_job = True;
        ## check if simulation have same parameters
        for s in list_vars_common:
            if not job_info[job1][s] == job_info[job2][s]:
                found_job = False;
        ## check if one simulation is parareal and the other is xbraid
        if p_or_x == "parareal":
            if not job_info[job2]["runtime.xbraid_enabled"]:
                found_job = False;
        elif p_or_x == "xbraid":
            if not job_info[job2]["runtime.parareal_enabled"]:
                found_job = False;

        if found_job:
            if p_or_x == "parareal":
                job_p = job1;
                job_x = job2;
            else:
                job_p = job2;
                job_x = job1;

            dt_p = job_info[job_p]["runtime.max_simulation_time"] / job_info[job_p]["runtime.parareal_coarse_slices"];
            dt_x = job_info[job_x]["runtime.timestep_size"] * job_info[job_x]["runtime.xbraid_cfactor"];

            if np.abs(dt_p - dt_x) > 1e-10:
                found_job = False;

        if found_job:
            ipair += 1;

            read_jobs.append(job2);

            if p_or_x == "parareal":
                list_files = glob.glob(path_simulations + "/" + job1 + "/parareal_error*");
            elif p_or_x == "xbraid":
                list_files = glob.glob(path_simulations + "/" + job1 + "/xbraid_error*");
            list_files = [os.path.basename(f) for f in list_files];

            max_diff = 0
            print("      -> Pair #{} : comparing {} files".format(ipair, len(list_files)));
            nb_not_found_files = 0;
            for f in list_files:

                if "_spec" not in f:

                    err_L1_1, err_L2_1, err_Linf_1 = read_error_file(path_simulations + "/" + job1 + "/" + f, error_type = "physical");
                    if p_or_x == "parareal":
                        fname2 = path_simulations + "/" + job2 + "/" + f.replace("parareal", "xbraid");
                    elif p_or_x == "xbraid":
                        fname2 = path_simulations + "/" + job2 + "/" + f.replace("xbraid", "parareal");
                    if not os.path.exists(fname2):
                        nb_not_found_files += 1;
                        continue;
                    err_L1_2, err_L2_2, err_Linf_2 = read_error_file(fname2, error_type = "physical");
                    err_Linf = np.abs(err_Linf_1 - err_Linf_2);
                    err_L1 = np.abs(err_L1_1 - err_L1_2);
                    err_L2 = np.abs(err_L2_1 - err_L2_2);
                    assert err_Linf < small, (err_Linf_1, err_Linf_2, np.abs(err_Linf_1 - err_Linf_2), f, job1, job2);
                    assert err_L1 < small, (err_L1_1, err_L1_2, f);
                    assert err_L2 < small, (err_L2_1, err_L2_2, f);
                    max_diff = np.max([max_diff, err_Linf, err_L1, err_L2]);
                else:
                    err_1 = read_error_file(path_simulations + "/" + job1 + "/" + f, error_type = "spectral");
                    if p_or_x == "parareal":
                        fname2 = path_simulations + "/" + job2 + "/" + f.replace("parareal", "xbraid");
                    elif p_or_x == "xbraid":
                        fname2 = path_simulations + "/" + job2 + "/" + f.replace("xbraid", "parareal");
                    if not os.path.exists(fname2):
                        nb_not_found_files += 1;
                        continue;
                    err_2 = read_error_file(fname2, error_type = "spectral");
                    for rnorm in err_1.keys():
                        err = np.abs(err_1[rnorm] - err_2[rnorm])
                        assert err < small, (rnorm, err_1, err_2, f)
                        max_diff = np.max([max_diff, err])

            print("     -> Number of files not found in both jobs: " + str(nb_not_found_files));
            print("     -> Max diff between errors: " + str(max_diff));
            print("                                             -> OK");
            break;
