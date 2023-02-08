#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule.JobMule import *
from itertools import product
from mule.utils import exec_program

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="exp_stability_plots"
jg.compile.mode = 'release'

jg.gen_jobscript_directory()


exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)


print("Benchmarks successfully finished")

print("Plotting results")

for i in [1, 2, 4]:

    # Errors
    file = "output_errors_plot_order_"+str(i)+".csv"
    exec_program(['mv', jg.job_dirpath+"/"+file, "./"], catch_output=False)
    print("Plotting "+file)
    exec_program(['./plot_stability_and_errors.py', file], catch_output=False)

    # Stability
    file = "output_stability_plot_order_"+str(i)+".csv"
    exec_program(['mv', jg.job_dirpath+"/"+file, "./"], catch_output=False)
    print("Plotting "+file)
    exec_program(['./plot_stability_and_errors.py', file], catch_output=False)


exec_program(['./tide_run_crop_and_convert_output_pdfs.sh'], catch_output=False)

#exec_program('mule.benchmark.cleanup_all', catch_output=False)

