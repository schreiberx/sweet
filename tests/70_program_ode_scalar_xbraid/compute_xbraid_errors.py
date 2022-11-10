#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob


def read_ref_solution(ref_path):

    ref_sol = {};

    list_ref_files = glob(ref_path + "/*csv");

    for f in list_ref_files:
        ## identify variable and time
        ff = os.path.basename(f).split("_t0");
        var = ff[0];
        t = float(ff[1].split(".csv")[0]);

        if not var in ref_sol.keys():
            ref_sol[var] = {};
        ref_sol[var][t] = np.loadtxt(f);

    return ref_sol;



list_vars = ["prog_h_pert", "prog_u", "prog_v"];
def read_xbraid_solution_compute_store_errors(path, ref_sol, ref_type):

    xbraid_sol = {};

    list_xbraid_files = glob(path + "/*csv");

    for f in list_xbraid_files:

        ## identify variable, time and iteration
        ff = os.path.basename(f).split("_t0");
        if ff[0][:8] == "residual":
            continue;
        var = ff[0];
        fff = ff[1].split("_iter");
        t = float(fff[0]);
        it = int(fff[1].split(".csv")[0]);

        ## check if csv file already contains computed errors
        if var[:14] == "xbraid_error":
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
        ####if t == 0:
        ####    continue;

        ref = ref_sol[var][t];
        sol = np.loadtxt(f);

        ###nx_ref, ny_ref = ref.shape;
        ###nx, ny = sol.shape;
        ###assert nx == nx_ref;
        ###assert ny == ny_ref;

        err_Abs = np.abs(sol - ref);
        err_Real = np.abs(sol.real - ref.real);
        err_Imag = np.abs(sol.imag - ref.imag);
        err_Phase = np.abs(np.angle(sol) - np.angle(ref));

        small = 1e-14;
        if np.abs(ref) > small:
            err_Abs /= np.abs(ref);
        if np.abs(ref.real) > small:
            err_Real /= np.abs(ref.real);
        if np.abs(ref.imag) > small:
            err_Imag /= np.abs(ref.imag);
        if np.abs(np.angle(ref)) > small:
            err_Phase /= np.abs(np.angle(ref));

        dirname = f.split("/")[0];

        error_file = open(dirname + "/xbraid_error_" + ref_type + "_" + os.path.basename(f)[7:], "w");
        error_file.write("errAbs {}\n".format(err_Abs));
        error_file.write("errReal {}\n".format(err_Real));
        error_file.write("errImag {}\n".format(err_Imag));
        error_file.write("errPhase {}".format(err_Phase));
        error_file.close();


##ref_sim = sys.argv[1];
fine_sim = sys.argv[1];

## get list of jobs in this directory
list_jobs = glob("job_bench_*");
## exclude fine simulation
list_jobs.remove(fine_sim);

## read ref and fine solutions
###ref_sol = read_ref_solution(ref_sim);
fine_sol = read_ref_solution(fine_sim);

## read xbraid solutions
for job in list_jobs:
    ###read_xbraid_solution_compute_store_errors(job, ref_sol, "ref");
    read_xbraid_solution_compute_store_errors(job, fine_sol, "fine");

