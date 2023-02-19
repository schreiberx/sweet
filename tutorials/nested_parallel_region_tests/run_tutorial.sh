#! /bin/bash




./benchmark_create_job_scripts.py

./compile_platform_*.sh

# This is just a compile benchmark because the OpenMP implementation might not work
# TO NOT ACTIVATE THIS PER DEFAULT!
#mule.benchmark.jobs_run_directly
