#! /usr/bin/env python3

import sys
import math

from mule.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Group together  similar time stepping methods
groups = ['runtime.xbraid_pt']

# Create plots for these variables
vars_ = ["phi_pert", "vrt", "div"]


###########################################################
# User configuration ends here ############################
###########################################################


tagnames_y = []
tag_cleanup_info = []

for i in vars_:
    tagnames_y += [
        f"sphere_data_diff_prog_{i}.res_norm_linf",
        f"sphere_data_diff_prog_{i}.res_norm_l1",
        f"sphere_data_diff_prog_{i}.res_norm_l2",
    ]

    tag_cleanup_info += [
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "res_norm_linf", "tag_dst": f"sphere_data_diff_prog_{i}.res_norm_linf"},
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "res_norm_l1", "tag_dst": f"sphere_data_diff_prog_{i}.res_norm_l1"},
        {"ref_file_starts_with": f"output_prog_{i}", "tag_src": "res_norm_l2", "tag_dst": f"sphere_data_diff_prog_{i}.res_norm_l2"},
    ]



j = JobsData(verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(" + "+key)

# Cleanup postprocessed data
JobsData_GroupsCleanupPostprocessed(job_groups, tag_cleanup_info, pickle_file_default_prefix="sphere_data_norms_physical_space_")

small = 1e-10

for tagname_y in tagnames_y:
    print("*"*80)
    print("Processing tagname "+tagname_y)
    print("*"*80)

    tagname_x = 'runtime.xbraid_pt'

    if True:
        """
        Use plotting format to create (x/y) data
        """
        d = JobsData_GroupsPlottingScattered(
                job_groups,
                tagname_x,
                tagname_y,
                meta_attribute_name = 'runtime.xbraid_nb_pt',
            )

        for group_name, group_data in d.get_data_float().items():

            print("*"*80)
            print("Nb. of processors: "+group_name)
            for (x, y) in zip(group_data['x_values'], group_data['y_values']):
                if abs(y) > small:
                    print("ERROR:", y)
                    raise Exception("Error should be zero!!!")
            print("OK!")

        print("*"*80)
        print("Test successful")
        print("*"*80)
