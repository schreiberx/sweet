#! /bin/bash

# Stop on first error
set -e

# First, cleanup things
mule.benchmark.cleanup_all

# Create all job directories
./1_benchmark_create_jobs.py

# Run benchmarks
mule.benchmark.jobs_run_directly

