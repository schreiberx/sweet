#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob
from read_bin_file import read_bin_file

file_type = "sweet"
list_vars = ["prog_phi_pert", "prog_vrt", "prog_div"];

def getArrayIndexByModes(n, m, N_max):

    assert n >= 0;
    assert n >= m;

    return (m * (2 * N_max - m + 1) >> 1)  + n;

def getMaxAbsRnorm(u, rnorm, N_max, verbose = False):

    err = 0;

    for m in range(int(rnorm)):
        for n in range(m, int(rnorm)):
            idx = getArrayIndexByModes(n, m, N_max);
            err = np.max([err, np.abs(u[idx] * np.conj(u[idx]))]);
            if verbose:
                print(m, n, idx, u[idx], err);

    return float(err);

def getMaxAbsRnormWithIndex(u, rnorm, N_max):

    err = 0;
    i = -1;

    for m in range(int(rnorm)):
        for n in range(m, int(rnorm)):
            idx = getArrayIndexByModes(n, m, N_max);
            err = np.max([err, np.abs(u[idx] * np.conj(u[idx]))]);
            if np.abs(err - np.abs(u[idx] * np.conj(u[idx]))) < 1e-13:
                i = idx;

    return float(err), i;



def read_ref_solution(ref_path):

    ref_sol = {};

    list_ref_files = glob(ref_path + "/*" + file_type);

    for f in list_ref_files:
        ## identify variable and time
        ff = os.path.basename(f).split("_t0");
        var = ff[0];
        t = float(ff[1].split("." + file_type)[0]);

        if not var in ref_sol.keys():
            ref_sol[var] = {};
        if file_type == "csv":
            ref_sol[var][t] = np.loadtxt(f)[1:, 1:];
        elif file_type == "sweet":
            ref_sol[var][t], m_max, n_max = read_bin_file(f);

    return ref_sol;



def read_parareal_solution_compute_store_errors(path, ref_sol, ref_type):

    parareal_sol = {};

    list_parareal_files = glob(path + "/*" + file_type);

    for f in list_parareal_files:

        ## identify variable, time and iteration
        ff = os.path.basename(f).split("_t0");
        if ff[0][:8] == "residual":
            continue;
        var = ff[0];
        fff = ff[1].split("_iter");
        t = float(fff[0]);
        it = int(fff[1].split("." + file_type)[0]);

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
        if file_type == "csv":
            sol = np.loadtxt(f)[1:, 1:];
        elif file_type == "sweet":
            sol, m_max, n_max = read_bin_file(f);

        if (file_type == "csv"):
            nx_ref, ny_ref = ref.shape;
            nx, ny = sol.shape;
            assert nx == nx_ref;
            assert ny == ny_ref;

            err_L1 = np.sum(np.abs(sol - ref)) / (nx * ny);
            err_L2 = np.sqrt(np.sum(np.abs(sol - ref)**2) / (nx * ny));
            err_Linf = np.max(np.abs(sol - ref));

            dirname = f.split("/")[0];

            error_file = open(dirname + "/parareal_error_" + ref_type + "_" + os.path.basename(f)[7:], "w");
            error_file.write("errL1 {}\n".format(err_L1));
            error_file.write("errL2 {}\n".format(err_L2));
            error_file.write("errLinf {}".format(err_Linf));
            error_file.close();

        elif (file_type == "sweet"):

            diff = sol - ref;

            n_modes = m_max + 1;

            rnorms = n_modes * np.array([1, 1./2., 1./4., 1./8., 1./16.]);

            dirname = f.split("/")[0];
            error_file = open(dirname + "/parareal_error_spec_" + ref_type + "_" + os.path.basename(f)[7:-5] + "csv", "w");
            ####print(dirname + "/parareal_error_spec_" + ref_type + "_" + os.path.basename(f)[7:-5] + "csv")

            eps_small = 1e-20;
            for rnorm in rnorms:
                norm_diff = getMaxAbsRnorm(diff, rnorm, n_modes - 1);
                norm_ref = getMaxAbsRnorm(ref, rnorm, n_modes - 1);
                if norm_diff < eps_small and norm_ref < eps_small:
                    err = 0.
                else:
                    err = norm_diff / norm_ref;
                error_file.write("errLinf {}".format(int(rnorm)) + " {}\n".format(err));
                ####err2, idx = getMaxAbsRnormWithIndex(diff, rnorm, n_modes - 1);
                ####print ("AAAA", var, it, t, rnorm, getMaxAbsRnorm(diff, rnorm, n_modes - 1), err2, idx, getMaxAbsRnorm(sol, rnorm, n_modes - 1), getMaxAbsRnorm(ref, rnorm, n_modes - 1), err);
                ####if "prog_div" in f and it == 4 and np.abs(t - 0.16) < 1e-10 and rnorm < 10:
                ####    print("DIFF");
                ####    getMaxAbsRnorm(diff, rnorm, n_modes - 1, verbose = True);
                ####    print("REF");
                ####    getMaxAbsRnorm(ref, rnorm, n_modes - 1, verbose = True);
                ####    print("SOL");
                ####    getMaxAbsRnorm(sol, rnorm, n_modes - 1, verbose = True);
            ####    print("DIFF");
            ####    for i in range(diff.size):
            ####        print (i, diff[i]);
            ####    print("REF");
            ####    for i in range(ref.size):
            ####        print (i, ref[i]);
            ####    print("SOL");
            ####    for i in range(sol.size):
            ####        print (i, sol[i]);

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

## read parareal solutions
for job in list_jobs:
    ###read_parareal_solution_compute_store_errors(job, ref_sol, "ref");
    read_parareal_solution_compute_store_errors(job, fine_sol, "fine");

