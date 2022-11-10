#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob

from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *


def read_error_file(path):

    lines = [line.rstrip() for line in open(path)];

    err_Abs = -1;
    err_Real = -1;
    err_Imag = -1;
    err_Phase = -1;

    for line in lines:
        spl = line.split();
        if spl[0] == "errAbs":
            err_Abs = float(spl[1]);
            continue;
        if spl[0] == "errReal":
            err_Real = float(spl[1]);
            continue;
        if spl[0] == "errImag":
            err_Imag = float(spl[1]);
            continue;
        if spl[0] == "errPhase":
            err_Phase = float(spl[1]);
            continue;

    if err_Abs < 0:
        print("ERROR: err_Abs not found in " + path);
        sys.exit();
    if err_Real < 0:
        print("ERROR: err_Real not found in " + path);
    if err_Imag < 0:
        print("ERROR: err_Imag not found in " + path);
    if err_Phase < 0:
        print("ERROR: err_Phase not found in " + path);

    return err_Abs, err_Real, err_Imag, err_Phase


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
small = 1e-10;
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

            list_files = glob.glob(path_simulations + "/" + job1 + "/xbraid_error*");
            list_files = [os.path.basename(f) for f in list_files];
            assert len(list_files) > 0

            max_diff = 0
            print("      -> Pair #{} : comparing {} files".format(ipair, len(list_files)));
            nb_not_found_files = 0;
            for f in list_files:
                err_Abs_1, err_Real_1, err_Imag_1, err_Phase_1 = read_error_file(path_simulations + "/" + job1 + "/" + f);
                if not os.path.exists(path_simulations + "/" + job2 + "/" + f):
                    nb_not_found_files += 1;
                    continue;
                err_Abs_2, err_Real_2, err_Imag_2, err_Phase_2 = read_error_file(path_simulations + "/" + job2 + "/" + f);
                ###print(err_L1_1, err_L1_2);
                ###print(err_L1_2, err_L2_2);
                ###print(err_Linf_1, err_Linf_2);
                ###print("");
                err_Abs = np.abs(err_Abs_1 - err_Abs_2);
                err_Real = np.abs(err_Real_1 - err_Real_2);
                err_Imag = np.abs(err_Imag_1 - err_Imag_2);
                err_Phase = np.abs(err_Phase_1 - err_Phase_2);
                assert err_Abs < small, (err_Abs_1, err_Abs_2, np.abs(err_Abs_1 - err_Abs_2), f, job1, job2);
                assert err_Real < small, (err_Real_1, err_Real_2, f);
                assert err_Imag < small, (err_Imag_1, err_Imag_2, f);
                assert err_Phase < small, (err_Phase_1, err_Phase_2, f);
                max_diff = np.max([max_diff, err_Abs, err_Real, err_Imag, err_Phase]);
            print("     -> Number of files not found in both jobs: " + str(nb_not_found_files));
            print("     -> Max diff between errors: " + str(max_diff));
            print("                                             -> OK");
            break;

