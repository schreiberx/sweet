#! /usr/bin/env python3

import numpy as np
import sys
import os
from glob import glob


def read_residuals(path):

    res = {};

    list_files = glob(path + "/residual_iter*.csv");

    print("Reading residual for ", len(list_files), "iterations");

    for f in list_files:
        print(f)
        ## identify variable and time
        ff = os.path.basename(f).split("_iter");
        niter = int(ff[1].split(".csv")[0]);

        fff = open(f)
        lines = fff.readlines();
        assert len(lines) == 4;
        spl = lines[3].split();
        assert len(spl) == 1;
        res[niter] = float(spl[0]);
        fff.close();

    return res;


def check_residuals(res, eps):

    res_max = 0;
    for niter in res.keys():
        assert(res[niter] < eps);
        res_max = max(res_max, res[niter]);

    print("Max residual: ", res_max);


## ugly check
## XBraid has no user function to get residuals at C-points, only print to terminal
## Then we get residuals from output.out
def read_check_residuals_C_points(path):
    f = open(path + "/output.out");

    lines = f.readlines();

    iline = 0;
    for line in lines:

        ## found line before printing C-point residuals
        if "Braid: || r_" in line:
            spl = line.split();
            niter = int(spl[2][2:]);

            ## check the follozing 2 * niter lines have residual equal to zero
            if len(lines[iline + 1]) > 1:
                print(" -- Checking C-points residuals at iter", niter);
                for j in range(iline + 1, iline + 1 + 2 * niter):
                    spl = lines[j].split();
                    res = float(spl[5]);
                    print(res)
                    assert res == 0.;
                print("    --> OK");

        iline += 1;

residual_type = sys.argv[1];
expected_max_residual = float(sys.argv[2]);

## get list of jobs in this directory
list_jobs = glob("job_bench_*");

if residual_type == "iteration":
    ## read parareal solutions
    for job in list_jobs:
        res = read_residuals(job);
        check_residuals(res, expected_max_residual);
elif residual_type == "C-point":
    for job in list_jobs:
        read_check_residuals_C_points(job);
else:
    sys.exit("Wrong residual type: " + residual_type);
