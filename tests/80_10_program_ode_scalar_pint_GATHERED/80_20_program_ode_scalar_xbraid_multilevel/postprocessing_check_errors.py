#! /usr/bin/env python3

import sys
import math

from mule.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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



# Group together  similar time stepping methods
groups = ['runtime.xbraid_max_levels', 'runtime.xbraid_cfactor', 'runtime.xbraid_pt']

# Create plots for these variables
vars_ = ["u"]


###########################################################
# User configuration ends here ############################
###########################################################


tagnames_y = []
tag_cleanup_info = []

for i in vars_:
    tagnames_y += [
        f"scalar_data_diff_prog_{i}.norm_l1",
        f"scalar_data_diff_prog_{i}.norm_l2",
        f"scalar_data_diff_prog_{i}.norm_linf"
    ]

    tag_cleanup_info += [
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "norm_linf", "tag_dst": f"scalar_data_diff_prog_{i}.norm_linf"},
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "norm_l1", "tag_dst": f"scalar_data_diff_prog_{i}.norm_l1"},
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "norm_l2", "tag_dst": f"scalar_data_diff_prog_{i}.norm_l2"}
    ]

j = JobsData(verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(" + "+key)


## check if there are two simulations per group (1 offline + 1 online)
for group in job_groups.values():
    assert len(group) == 2

## Cleanup postprocessed data
JobsData_GroupsCleanupPostprocessed(job_groups, tag_cleanup_info, pickle_file_default_prefix="scalar_data_norms_physical_space_", pint = True)

small = 1e-15
## compare online and offline errors
for key, group in job_groups.items():

    print(" ~ Group: ", key)

    niter = 0;
    while True:

        print('--> Comparing solutions at iteration ' + str(niter))

        niter_str = "_iter" + str(niter).zfill(3)

        for i in vars_:

            print(' == Var: ' + i)

            err_online = {}
            err_offline = {}
            ## identify offline and online
            for job in group.keys():
                if group[job]['runtime.xbraid_store_iterations']:
                    offline_job = job
                else:
                    online_job = job

            br = False

            ## get offline errors
            try:
                err_offline['norm_linf'] = group[offline_job][f"scalar_data_diff_prog_{i}.norm_linf" + niter_str]
                err_offline['norm_l1'] = group[offline_job][f"scalar_data_diff_prog_{i}.norm_l1" + niter_str]
                err_offline['norm_l2'] = group[offline_job][f"scalar_data_diff_prog_{i}.norm_l2" + niter_str]
            except:
                br = True
                break

            ## identify files to read online error
            error_files = mule.postprocessing.utils.get_job_output_files(group[offline_job])
            error_files = [f.replace("output_", "xbraid_error_fine_") for f in error_files if niter_str in f and f"prog_{i}" in f]

            path = group[online_job]['runtime.p_job_dirpath']
            for f in error_files:
                err_online['norm_l1'], err_online['norm_l2'], err_online['norm_linf'] = read_error_file(path + "/" + f)

            max_error = 0
            max_norm = 0
            for n in ['norm_l1', 'norm_l2', 'norm_linf']:
                assert abs(err_online[n] - err_offline[n]) < small, (err_online[n], err_offline[n])
                max_error = max(max_error, abs(err_online[n] - err_offline[n]))
                max_norm = max(max_norm, abs(err_online[n]), abs(err_offline[n]))
            print("   * Max norm: " + str(max_norm))
            print("   * Max error: " + str(max_error) + " (OK!)")

        if br:
            break
        niter += 1;
        print("")

    print("\n\n")

print("*"*80)
print("Test successful")
print("*"*80)
