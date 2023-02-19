#! /bin/bash

# Stop on first error
set -e

# Create all job directories
./benchmark_create_jobs.py

# Run benchmarks
mule.benchmark.jobs_run_directly_nonstop

# Compute all norms of difference between job output data and reference job(s)
mule.postprocessing.pickle.alljobs.sphere_data_norms_physical_space

# Plot nice looking pictures
./postprocessing_plot.py
