#!/bin/bash

get_tsm(){
	tsm=$1
	tsm=${tsm#*TS_};
	tsm=${tsm%.hpp};
	echo "$tsm";
}

dirname="simulations_offline_error";

echo_info "Cleaning up..."
mule.benchmark.cleanup_all || exit 1
rm -r $dirname;
mkdir simulations_offline_error;

for online_error in {0,1};do
	for tsm_fine in ../../src/programs/swe_plane_timeintegrators/SWE_Plane_TS*hpp; do

		tsm_fine=$(get_tsm $tsm_fine);
		if [ "$tsm_fine" = "interface" ]; then
		  continue;
		fi

		for tsm_coarse in ../../src/programs/swe_plane_timeintegrators/SWE_Plane_TS*hpp; do

			tsm_coarse=$(get_tsm $tsm_coarse);
			if [ "$tsm_coarse" = "interface" ]; then
			  continue;
			fi

			echo $tsm_fine $tsm_coarse
			dirname2=simulations_offline_error"/"${tsm_fine}"_"${tsm_coarse}

			## parareal tests without online error computation
			echo_info "Creating parareal simulations with:"
			echo_info "tsm_fine" $tsm_fine;
			echo_info "tsm_coarse" $tsm_coarse;
			./benchmarks_create.py $tsm_fine $tsm_coarse $online_error $ref_sim $dirname2"/"$fine_sim > dummy || exit 1

			## identify ref simulation
			ref_sim=$(cat ref_sim);
			echo REF simulation: $ref_sim;

			## identify fine simulation
			fine_sim=$(cat fine_sim);
			echo FINE simulation : $fine_sim;

			echo_info "Running simulations"
			mule.benchmark.jobs_run_directly > dummy|| exit 1

			echo_info "Computing errors"
			./compute_parareal_errors.py $ref_sim $fine_sim || exit 1

			if [ $online_error -eq 0 ]; then
				echo_info "Moving simulations into " $dirname2;
				mkdir $dirname2;
				mv job_bench_* $dirname2;
			fi;

			rm ref_sim;
			rm fine_sim;

		done;
	done;
done;
