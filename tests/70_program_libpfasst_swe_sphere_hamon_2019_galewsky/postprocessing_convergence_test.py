#! /usr/bin/env python3

import sys
import math

from mule_local.JobMule import *
from mule.plotting.Plotting import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



groups = [
    'runtime.libpfasst_nlevels',
    'runtime.libpfasst_niters',
    'runtime.libpfasst_nnodes',
    'runtime.libpfasst_coarsening_multiplier'
]


tagnames_y = [
    'sphere_data_diff_prog_vrt.res_norm_linf',
]


j = JobsData(verbosity=0)

c = JobsDataConsolidate(j)
print("")
print("Groups:")
job_groups = c.create_groups(groups)
for key, g in job_groups.items():
    print(" + "+key)


for tagname_y in tagnames_y:
    print("*"*80)
    print("Processing tagname "+tagname_y)
    print("*"*80)

    tagname_x = 'runtime.timestep_size'

    if True:
        """
        Use plotting format to create (x/y) data
        """
        d = JobsData_GroupsPlottingScattered(
                job_groups,
                tagname_x,
                tagname_y,
                meta_attribute_name = 'runtime.timestepping_order',
            )

        for group_name, group_data in d.get_data_float().items():
            print("*"*80)
            print("Group: "+group_name)
            prev_value = -1.0
            conv = '-'
            convergence_order = None
            for (x, y, convergence_order_) in zip(group_data['x_values'], group_data['y_values'], group_data['meta_values']):

                if prev_value > 0:
                    conv = y/prev_value
                elif prev_value == 0:
                    conv = '[error=0]'

                print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv))
                prev_value = y

                if convergence_order == None:
                    convergence_order = convergence_order_
                else:
                    if convergence_order != convergence_order_:
                        raise Exception("Convergence order mismatch!!!")


            print("")
            print("Testing convergence")

            #
            # Setup default values
            #

            # 'convergence', 'error'
            test_type = 'convergence'

            # Convergence tolerance
            error_tolerance_convergence = 0.15

            # Range to check for convergence
            conv_test_range_start = 0
            conv_test_range_end = 4

            if 'vrt' in tagname_y or 'div' in tagname_y:

                if 'exp' in group_name:
                    test_type = 'error'
                    error_tolerance_error = 1e-12


            elif 'phi' in tagname_y:
                if 'exp' in group_name:
                    test_type = 'error'
                    error_tolerance_error = 1e-4

            else:
                raise Exception("Tagname "+tagname_y+" unknown")


            print(" + test_type: "+test_type)

            if test_type == 'convergence':
                print(" + error_tolerance_convergence: "+str(error_tolerance_convergence))
            elif test_type == 'error':
                print(" + error_tolerance_error: "+str(error_tolerance_error))

            print(" + range start/end: "+str(conv_test_range_start)+", "+str(conv_test_range_end))

            if len(group_data['meta_values']) < conv_test_range_end:
                raise Exception("Not enough samples to run convergence test")

            for i in range(len(group_data['meta_values'])):
                if group_data['meta_values'][i] != group_data['meta_values'][0]:
                    print("FATAL: Different convergence orders in same test")
                    for i in range(len(group_data['meta_values'])):
                        print("order: "+str(group_data['meta_values']))

                    raise Exception("FATAL: Different convergence orders in same test")

            l = len(group_data['x_values'])
            if l < conv_test_range_end:
                print("There are only "+str(l)+" values, but we need at least "+str(conv_test_range_end)+" values")
                raise Exception("Not enough values to study convergence")

            prev_value = -1.0
            conv = '-'
            for i in range(conv_test_range_start, conv_test_range_end):
                x = group_data['x_values'][i]
                y = group_data['y_values'][i]
                meta = group_data['meta_values'][i]

                if prev_value > 0:
                    conv = y/prev_value
                elif prev_value == 0:
                    conv = '[error=0]'

                error_convergence = '-'
                if isinstance(conv, float):
                    # Convergence order is stored in meta value
                    target_conv = pow(2.0, meta)
                    error_convergence = abs(conv - target_conv)/target_conv

                print("\t"+str(x)+"\t=>\t"+str(y)+"\tconvergence: "+str(conv)+"\terror: "+str(error_convergence))

                if test_type == 'convergence':
                    # Test for convergence if exists
                    if error_convergence != '-':
                        if error_convergence > error_tolerance_convergence:
                            print("Error: "+str(error_convergence))
                            if len(sys.argv) <= 1:
                                raise Exception("Convergence exceeds tolerance of "+str(error_tolerance_convergence))

                elif test_type == 'error':
                    # Alternate tests instead of convergence check
                    # Convergence doesn't really make sense for REXI in the way how it's applied
                    # This should be only used for l_exp and lg_exp
                    # Just ensure that the errors are below a certain level
                    if y > error_tolerance_error:
                        print("Error: "+str(y))
                        if len(sys.argv) <= 1:
                             raise Exception("Error exceeds tolerance of "+str(error_tolerance_error))

                else:
                    raise Exception("Unknown test type "+test_type)

                prev_value = y

            if len(sys.argv) <= 1:
                print("[OK]")

        if len(sys.argv) <= 1:
            print("*"*80)
            print("Convergence tests successful")
            print("*"*80)
