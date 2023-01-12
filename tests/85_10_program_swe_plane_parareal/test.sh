#!/bin/bash

get_tsm(){
	tsm=$1
	tsm=${tsm#*TS_};
	tsm=${tsm%.hpp};
	echo "$tsm";
}

set -e

dirname="simulations_offline_error";

cd "$(dirname $0)"

echo_info "Cleaning up..."
mule.benchmark.cleanup_all || exit 1
if [ -d $dirname ]; then
	rm -r $dirname;
fi
mkdir $dirname;


echo ""

for i in {0,1,2};do
	##for tsm_fine in ../../src/programs/swe_plane_timeintegrators/SWE_Plane_TS*hpp; do ## full version
	for tsm_fine in {l_cn,l_erk,l_cn_na_sl_nd_settls,l_rexi_n_erk,l_rexi_na_sl_nd_etdrk}; do ## short version
	##for tsm_fine in l_cn; do ## short version

		tsm_fine=$(get_tsm $tsm_fine);
		if [ "$tsm_fine" = "interface" ]; then
		  continue;
		fi

		##for tsm_coarse in ../../src/programs/swe_plane_timeintegrators/SWE_Plane_TS*hpp; do ## full version
		for tsm_coarse in {l_cn,l_erk,l_cn_na_sl_nd_settls,l_rexi_n_erk,l_rexi_na_sl_nd_etdrk}; do ## short version
		##for tsm_coarse in l_cn; do ## short version

			tsm_coarse=$(get_tsm $tsm_coarse);
			if [ "$tsm_coarse" = "interface" ]; then
			  continue;
			fi

			dirname2=simulations_offline_error"/"${tsm_fine}"_"${tsm_coarse}

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
				ref_sim=$(cat ref_sim);

				## identify fine simulation
				fine_sim=$(cat fine_sim);

				mv $dirname2/job_bench* .;

				echo_info "---> Computing errors with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse
				./compute_parareal_errors.py $ref_sim $fine_sim || exit 1

                                mv ref_sim $dirname2/.;
                                mv fine_sim $dirname2/.;
                        fi;

			## only parareal with online error computation
			if [ $i == 2 ]; then
				## parareal tests without online error computation
				echo_info "---> Running parareal simulations (online error computation) with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

				## identify ref simulation
				ref_sim=$(cat $dirname2/ref_sim);

				## identify fine simulation
				fine_sim=$(cat $dirname2/fine_sim);


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

		done
	done
done

mule.benchmark.cleanup_all || exit 1

rm -r "$dirname"
rm -r "$dirname2"
rm -f tmp_job_benchmark_create_dummy.txt

echo ""
echo_info "Test successful!"
