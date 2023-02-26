#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob


def read_ref_solution(ref_path):

    ref_sol = {};

    list_ref_files = glob(ref_path + "/*csv");

    for f in list_ref_files:
        if "amp_phase" in f:
            continue;
        ## identify variable and time
        ff = os.path.basename(f).split("_t0");
        var = ff[0];
        t = float(ff[1].split(".csv")[0]);

        if not var in ref_sol.keys():
            ref_sol[var] = {};
        ref_sol[var][t] = np.loadtxt(f);

    return ref_sol;



list_vars = ["prog_u", "prog_v"];
def read_parareal_solution_compute_store_errors(path, ref_sol, ref_type):

    parareal_sol = {};

    list_parareal_files = glob(path + "/*csv");

    ###print("LIST FILES", list_parareal_files)
    for f in list_parareal_files:
        if "amp_phase" in f:
            continue;
        if "_spec_" in f:
            continue;

        ## identify variable, time and iteration
        ff = os.path.basename(f).split("_t0");
        var = ff[0];
        if "_iter" in ff:
            fff = ff[1].split("_iter");
            t = float(fff[0]);
            it = int(fff[1].split(".csv")[0]);
        else:
            foo = ff[1].split(".csv")
            t = float(foo[0]);
            it = -1

        #print(f"t: {t}")
        #print(f"it: {it}")

        ## check if csv file already contains computed errors
        if var[:14] == "parareal_error":
            continue;

        ## check if csv file contains listed vars
        var_ok = False;
        for v in list_vars:
            if v in var:
                var_ok = True;
                break;
        if not var_ok:
            continue;

        ## check if not initial solution
        if t == 0:
            continue;

        ref = ref_sol[var][t];
        sol = np.loadtxt(f);

        nx_ref, ny_ref = ref.shape;
        nx, ny = sol.shape;
        assert nx == nx_ref;
        assert ny == ny_ref;

        err_L1 = np.sum(np.abs(sol - ref)) / (nx * ny);
        err_L2 = np.sqrt(np.sum(np.abs(sol - ref)**2) / (nx * ny));
        err_Linf = np.max(np.abs(sol - ref));

        dirname = f.split("/")[0];

        errorfilepath = dirname + "/parareal_error_" + ref_type + "_" + os.path.basename(f)[7:]

        print(f"Writing to error file '{errorfilepath}'")
        error_file = open(errorfilepath, "w");
        error_file.write("errL1 {}\n".format(err_L1));
        error_file.write("errL2 {}\n".format(err_L2));
        error_file.write("errLinf {}".format(err_Linf));
        error_file.close();


ref_sim = sys.argv[1];
fine_sim = sys.argv[2];

## get list of jobs in this directory
list_jobs = glob("job_bench_*");
## exclude fine simulation
list_jobs.remove(fine_sim);

## read ref and fine solutions
ref_sol = read_ref_solution(ref_sim);
fine_sol = read_ref_solution(fine_sim);

## read parareal solutions
for job in list_jobs:
    read_parareal_solution_compute_store_errors(job, ref_sol, "ref");
    read_parareal_solution_compute_store_errors(job, fine_sol, "fine");

