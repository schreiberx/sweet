#!/bin/bash

cd "$(dirname $0)"

get_tsm(){
	tsm=$1
	tsm=${tsm#*TS_};
	tsm=${tsm%.hpp};
	echo "$tsm";
}

set -e


dirname_="output_simulations_offline_error";

echo_info "Cleaning up..."
mule.benchmark.cleanup_all || exit 1
if [ -d $dirname_ ]; then
	rm -rf $dirname_;
fi
mkdir $dirname_;


echo ""

for i in {0,1,2};do
	for tsm_fine in dummy; do ## short version

		tsm_fine=$(get_tsm $tsm_fine);
		if [ "$tsm_fine" = "interface" ]; then
		  continue;
		fi

		for tsm_coarse in dummy; do ## short version

			tsm_coarse=$(get_tsm $tsm_coarse);
			if [ "$tsm_coarse" = "interface" ]; then
			  continue;
			fi

			dirname2=output_simulations_offline_error"/"${tsm_fine}"_"${tsm_coarse}

			## only parareal
			if [ $i == 0 ]; then
				## parareal tests without online error computation
				echo_info "---> Running parareal simulations (offline error computation) with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

				./benchmarks_create.py $tsm_fine $tsm_coarse parareal 0 $ref_sim $dirname2"/"$fine_sim > tmp_job_benchmark_create_dummy.txt || exit 1

				mule.benchmark.jobs_run_directly || exit 1
                        fi;

			## fine and ref
			if [ $i == 1 ]; then
				## parareal tests without online error computation
				echo_info "---> Running fine and ref simulations with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

				./benchmarks_create.py $tsm_fine $tsm_coarse ref 0 $ref_sim $dirname2"/"$fine_sim  > tmp_job_benchmark_create_dummy.txt || exit 1

				mule.benchmark.jobs_run_directly|| exit 1

				## identify ref simulation
				ref_sim=$(cat tmp_ref_sim.txt);

				## identify fine simulation
				fine_sim=$(cat tmp_fine_sim.txt);

				mv $dirname2/job_bench* .;

				echo_info "---> Computing errors with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse
				./compute_parareal_errors.py $ref_sim $fine_sim || exit 1

                                mv tmp_ref_sim.txt $dirname2/.;
                                mv tmp_fine_sim.txt $dirname2/.;
                        fi;

			## only parareal with online error computation
			if [ $i == 2 ]; then
				## parareal tests without online error computation
				echo_info "---> Running parareal simulations (online error computation) with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

				## identify ref simulation
				ref_sim=$(cat $dirname2/tmp_ref_sim.txt)

				## identify fine simulation
				fine_sim=$(cat $dirname2/tmp_fine_sim.txt)

				./benchmarks_create.py $tsm_fine $tsm_coarse parareal 1 ../$dirname2"/"$ref_sim ../$dirname2"/"$fine_sim > tmp_job_benchmark_create_dummy.txt || exit 1

				mv $ref_sim $dirname2/.

				mule.benchmark.jobs_run_directly || exit 1

				mv $dirname2/$ref_sim .

                        fi;

			if [ $i -eq 0 ]; then
				mkdir $dirname2;
			fi;
			if [ $i -le 2 ]; then
				mv job_bench_* $dirname2;
			fi;
			if [ $i -eq 2 ]; then
				echo_info "---> Comparing online and offline errors with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse
				./compare_online_offline_errors.py $dirname2 $fine_sim
			fi;
			echo ""

		done;
	done;
done;

rm -rf $dirname_
rm -rf $dirname2
rm -f tmp_job_benchmark_create_dummy.txt

mule.benchmark.cleanup_all || exit 1

echo ""
echo_info "Test successful!"
