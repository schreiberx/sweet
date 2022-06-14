#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob

from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *


def read_error_file(path):

    lines = [line.rstrip() for line in open(path)];

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
        print("ERROR: err_L1 not found in " + path);
        sys.exit();
    if err_L2 < 0:
        print("ERROR: err_L2 not found in " + path);
    if err_Linf < 0:
        print("ERROR: err_Linf not found in " + path);

    return err_L1, err_L2, err_Linf


path_simulations = sys.argv[1];
fine_sim = sys.argv[2];

## get list of jobs in this directory
list_jobs = glob.glob(path_simulations + "/job_bench_*");
list_jobs = [os.path.basename(job) for job in list_jobs];
## exclude fine simulation
list_jobs.remove(fine_sim);

print("    ** {} jobs found.".format(len(list_jobs)));

jd = JobsData(path_simulations + '/job_bench_*', verbosity=0).get_flattened_data();
## get useful job info
job_info = {};
for key in jd.keys():
    path = os.path.basename(jd[key]["jobgeneration.p_job_dirpath"]);
    if path == fine_sim:
        continue;
    job_info[path] = {};
    for s in ["runtime.xbraid_cfactor", "runtime.xbraid_max_levels", "runtime.xbraid_pt", "runtime.xbraid_store_iterations"]:
        job_info[path][s] = jd[key][s];


## find identical jobs
small = 1e-8;
read_jobs = [];
ipair = 0;
for job1 in list_jobs:
    if job1 in read_jobs:
        continue;
    read_jobs.append(job1);

    found_job = False;
    for job2 in list_jobs:
        if job2 in read_jobs:
            continue;
        found_job = True;
        for s in ["runtime.xbraid_cfactor", "runtime.xbraid_max_levels", "runtime.xbraid_pt"]:
            if not job_info[job1][s] == job_info[job2][s]:
                found_job = False;


        if found_job:
            ipair += 1;
            assert not job_info[job1]["runtime.xbraid_store_iterations"] == job_info[job2]["runtime.xbraid_store_iterations"];

            read_jobs.append(job2);

            list_files = glob.glob(path_simulations + "/" + job1 + "/parareal_error*");
            list_files = [os.path.basename(f) for f in list_files];

            max_diff = 0
            print("      -> Pair #{} : comparing {} files".format(ipair, len(list_files)));
            for f in list_files:
                err_L1_1, err_L2_1, err_Linf_1 = read_error_file(path_simulations + "/" + job1 + "/" + f);
                err_L1_2, err_L2_2, err_Linf_2 = read_error_file(path_simulations + "/" + job2 + "/" + f);
                ###print(err_L1_1, err_L1_2);
                ###print(err_L1_2, err_L2_2);
                ###print(err_Linf_1, err_Linf_2);
                ###print("");
                err_Linf = np.abs(err_Linf_1 - err_Linf_2);
                err_L1 = np.abs(err_L1_1 - err_L1_2);
                err_L2 = np.abs(err_L2_1 - err_L2_2);
                assert err_Linf < small, (err_Linf_1, err_Linf_2, np.abs(err_Linf_1 - err_Linf_2), f, job1, job2);
                assert err_L1 < small, (err_L1_1, err_L1_2, f);
                assert err_L2 < small, (err_L2_1, err_L2_2, f);
                max_diff = np.max([max_diff, err_Linf, err_L1, err_L2]);
            print("     -> Max diff between errors: " + str(max_diff));
            print("                                             -> OK");
            break;

