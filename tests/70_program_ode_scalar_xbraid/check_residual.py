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


expected_max_residual = float(sys.argv[1]);

## get list of jobs in this directory
list_jobs = glob("job_bench_*");

## read parareal solutions
for job in list_jobs:
    res = read_residuals(job);
    check_residuals(res, expected_max_residual);

