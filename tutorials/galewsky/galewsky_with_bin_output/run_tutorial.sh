#! /bin/bash

set -e

./0_clean.sh
./1_benchmark_create_jobs.py
./2_benchmark_compile.sh
./3_benchmark_run_all.sh
./4_postprocess_plots.sh

