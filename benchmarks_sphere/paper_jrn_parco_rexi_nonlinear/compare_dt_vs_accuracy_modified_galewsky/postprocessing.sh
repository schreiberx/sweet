#! /usr/bin/env bash



#
# Step 1)
# Compute L1,2,inf and RMS norms
# and store them to
#	[job]/sphere_data_diff.pickle
#


if false; then
#if true; then
	REF_JOB="$(ls -1 -d ./job_benchref_*/)"
	CMP_JOBS="$(ls -1 -d ./job_bench_*/)"

	if [[ -z "$REF_JOB" ]]; then
		echo_error_exit "No reference job found"
	fi

	./pp_compute_norms_refjob_to_jobs.sh 		\
		$REF_JOB				\
		t00000000120.00000000.csv		\
		$CMP_JOBS
fi

#if false; then
if true; then
	./pp_plot_dt_vs_accuracy.py sphere_data_diff_prog_h.norm_l1
fi

