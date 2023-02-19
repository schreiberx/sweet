#! /bin/bash

# Stop on first error
set -e

# First, cleanup things
mule.benchmark.cleanup_all

# Create all job directories
./1_benchmark_create_jobs.py

# Run benchmarks
mule.benchmark.jobs_run_directly

# Compute all norms of difference between job output data and reference job(s)
mule.postprocessing.pickle.alljobs.sphere_data_norms_physical_space

# Plot nice looking pictures
./2_postprocessing_plot.py
