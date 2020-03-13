#! /usr/bin/env bash


mule.benchmark.cleanup_all

./benchmarks_create.py || exit 1

mule.benchmark.jobs_run_directly || exit 1

./postprocessing.py || exit 1

mule.benchmark.cleanup_all
