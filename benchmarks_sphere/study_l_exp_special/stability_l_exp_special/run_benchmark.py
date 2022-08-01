#! /usr/bin/env python3

import sys
import os
os.chdir(os.path.dirname(sys.argv[0]))

from mule_local.JobMule import *
from itertools import product
from mule.exec_program import *

exec_program('mule.benchmark.cleanup_all', catch_output=False)

jg = JobGeneration()
jg.compile.program="exp_stability_plots"
jg.compile.mode = 'release'

jg.gen_jobscript_directory()

exitcode = exec_program('mule.benchmark.jobs_run_directly', catch_output=False)
if exitcode != 0:
    sys.exit(exitcode)

print("Benchmarks successfully finished")

exec_program('mule.benchmark.cleanup_all', catch_output=False)


for i in [1, 2, 4]:
    file = "output_stability_plot_order"+str(i)+".csv"
    print("Plotting "+file)

    exec_program(['./plot_stability.py', file], catch_output=False)

