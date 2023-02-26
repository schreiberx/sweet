#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## itest = 0: test user-defined wrapped functions
## itest = 1: max_levels = 1 -> solution should be equal to the serial one
## itest = 2: max_levels = 1 and 2 processors in time -> idem
## itest = 3: max_levels = 2, tol = 0., max_iter = 3, use_seqsoln, n processors in time -> residul should be zero at each iteration
## itest = 4: print_level = 3 -> check residual norm at each C-point: should be zero for two C-points per iteration
## itest = 5: multilevel tests + n processors in time -> XBraid solution within tol w.r.t. serial solution
###############
## Additional unit tests:
## itest = 6: basic test for spatial coarsening feature: simulations with spatial_coarsening = 0 or = N ( = fine resolution) should be identical
## itest = 7: check if parareal and xbraid with specific configurations provide identical results
###############


get_tsm(){
	tsm=$1
	tsm=${tsm#*TS_}
	tsm=${tsm%.hpp}
	echo "$tsm"
}

set -e

dirname_serial="tmp_serial"
dirname_offline_error="output_simulations_offline_error"

cd "$(dirname $0)"

echo_info "Cleaning up..."
mule.benchmark.cleanup_all || exit 1



echo ""

for itest in {-1..7};do
	echo "*********************"
	echo "Running debug test" $itest
	echo "*********************"

	echo "Description:"
	if [ "$itest" == -1  ]; then
		echo "Serial simulation"
	elif [ "$itest" == 0  ]; then
		echo "Wrapper tests"
	elif [ "$itest" == 1 ]; then
		echo "Max levels = 1"
	elif [ "$itest" == 2 ]; then
		echo "Max levels = 1 + 2 processors in time"
	elif [ "$itest" == 3 ]; then
		echo "Max levels = 2, tol = 0, max_iter = 3, use_seqsoln"
	elif [ "$itest" == 4 ]; then
		echo "print_level = 3"
	elif [ "$itest" == 5 ]; then
		echo "Misc. multilevel tests"
	fi
	echo ""

	tsm_fine="ln_erk"
	tsm_coarse="ln_erk"


	if [ "$itest" == -1  ]; then
		mkdir "$dirname_serial"

		./benchmarks_create.py ref $itest $tsm_fine $tsm_coarse 1 > tmp_job_benchmark_create_dummy.txt || exit 1
		time mule.benchmark.jobs_run_directly || exit 1
		mv job_bench_* "$dirname_serial/"

		# Get job directory name for this reference solution
		# This will be reused throughout all other test cases
		fine_sim1=$(cat tmp_fine_sim.txt)

		mule.benchmark.cleanup_job_dirs || exit 1

	elif [ "$itest" == 0 ]; then
		./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse 1 > tmp_job_benchmark_create_dummy.txt || exit 1
		time mule.benchmark.jobs_run_directly || exit 1
		mule.benchmark.cleanup_job_dirs || exit 1

	elif [ "$itest" == 1 ] || [ "$itest" == 2 ]; then
		./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $itest > tmp_job_benchmark_create_dummy.txt || exit 1
		time mule.benchmark.jobs_run_directly || exit 1
		cp -r "$dirname_serial/$fine_sim1" .
		./compare_to_fine_solution.py $fine_sim1
		mule.benchmark.cleanup_job_dirs || exit 1

	elif [ "$itest" == 3 ]; then
		for nproc in {1,2}; do
			echo "  -------------"
			echo "  -- nproc:" $nproc
			echo "  -------------"
			./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $nproc > tmp_job_benchmark_create_dummy.txt || exit 1
			time mule.benchmark.jobs_run_directly || exit 1
			./check_residual.py iteration 1e-16
			mule.benchmark.cleanup_job_dirs || exit 1
			echo ""
		done

	elif [ "$itest" == 4 ]; then
		for nproc in {1,2}; do
			echo "  -------------"
			echo "  -- nproc:" $nproc
			echo "  -------------"
			./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $nproc > tmp_job_benchmark_create_dummy.txt || exit 1
			time mule.benchmark.jobs_run_directly || exit 1
			./check_residual.py C-point 1e-16
			mule.benchmark.cleanup_job_dirs || exit 1
			echo ""
		done
	elif [ "$itest" == 5 ]; then

		for nproc in {1,2}; do

			echo "  -------------"
			echo "  -- nproc:" $nproc
			echo "  -------------"

			mule.benchmark.cleanup_job_dirs || exit 1


			##for tsm_fine in ../../src/programs/swe_plane_timeintegrators/SWE_Plane_TS*hpp; do ## full version
			#for tsm_fine in {l_cn,l_erk,l_cn_na_sl_nd_settls,l_rexi_n_etdrk,l_rexi_na_sl_nd_etdrk}; do ## short version
			for tsm_fine in {l_cn,l_erk,l_cn_na_sl_nd_settls}; do ## other short version
			#for tsm_fine in {l_cn,l_rexi_n_etdrk}; do ## short version

				tsm_fine=$(get_tsm $tsm_fine)
				if [ "$tsm_fine" = "interface" ]; then
					continue
				fi

				##for tsm_coarse in ../../src/programs/swe_plane_timeintegrators/SWE_Plane_TS*hpp; do ## full version
				#for tsm_coarse in {l_cn,l_erk,l_cn_na_sl_nd_settls,l_rexi_n_etdrk,l_rexi_na_sl_nd_etdrk}; do ## short version
				for tsm_coarse in {l_cn,l_erk,l_cn_na_sl_nd_settls}; do ## other short version
				#for tsm_coarse in {l_cn,l_rexi_n_etdrk}; do ## short version

					tsm_coarse=$(get_tsm $tsm_coarse)
					if [ "$tsm_coarse" = "interface" ]; then
						continue
					fi

					mkdir "$dirname_offline_error"

					for i in {0,1,2}; do

						dirname2=${dirname_offline_error}"/"${tsm_fine}"_"${tsm_coarse}

						## only xbraid
						if [ $i == 0 ]; then
							## xbraid tests without online error computation
							echo_info "---> Running XBraid simulations (offline error computation) with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

							./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $nproc 0 $dirname2"/"$fine_sim1 > tmp_job_benchmark_create_dummy.txt || exit 1

							time mule.benchmark.jobs_run_directly || exit 1
						fi

						## fine and ref
						if [ $i == 1 ]; then
							## parareal tests without online error computation
							echo_info "---> Running fine and ref simulations with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

							./benchmarks_create.py ref $itest $tsm_fine $tsm_coarse 1 0 $dirname2"/"$fine_sim1  > tmp_job_benchmark_create_dummy.txt || exit 1

							time mule.benchmark.jobs_run_directly|| exit 1

							## identify fine simulation
							fine_sim2=$(cat tmp_fine_sim.txt)

							mv $dirname2/job_bench* .

							echo_info "---> Computing errors with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse
							./compute_xbraid_errors.py $fine_sim2 || exit 1

							mv tmp_fine_sim.txt $dirname2/.

						fi

						## only xbraid with online error computation
						if [ $i == 2 ]; then
							echo_info "---> Running XBraid simulations (online error computation) with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse

							## identify fine simulation
							fine_sim2=$(cat $dirname2/tmp_fine_sim.txt)

							./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $nproc 1 ../$dirname2"/"$fine_sim2 > tmp_job_benchmark_create_dummy.txt || exit 1

							time mule.benchmark.jobs_run_directly || exit 1

						fi

						if [ $i -eq 0 ]; then
							mkdir "$dirname2"
						fi
						if [ $i -le 2 ]; then
							mv job_bench_* "$dirname2"
						fi
						if [ $i -eq 2 ]; then
							echo_info "---> Comparing online and offline errors with tsm_fine and tsm_coarse:" $tsm_fine $tsm_coarse
							./compare_online_offline_errors.py $dirname2 $fine_sim2
						fi
						echo ""

					done
			
					mule.benchmark.cleanup_job_dirs || exit 1
					rm -r "$dirname_offline_error"

				done # coarse
			done # fine


		done # Nproc

	elif [ "$itest" == 6 ]; then

		for nproc in {1,2}; do

			echo "  -------------"
			echo "  -- nproc:" $nproc
			echo "  -------------"

			./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $nproc 1 ../$dirname_serial"/"$fine_sim1 > tmp_job_benchmark_create_dummy.txt || exit 1

			time mule.benchmark.jobs_run_directly || exit 1
			./compare_online_offline_errors.py . $fine_sim1 1

			mule.benchmark.cleanup_job_dirs || exit 1
		done

	elif [ "$itest" == 7 ]; then

		for nproc in {1,5}; do

			echo "  -------------"
			echo "  -- nproc:" $nproc
			echo "  -------------"

			./benchmarks_create.py xbraid $itest $tsm_fine $tsm_coarse $nproc 1 ../$dirname_serial"/"$fine_sim1 > tmp_job_benchmark_create_dummy.txt || exit 1

			time mule.benchmark.jobs_run_directly || exit 1
			./compare_parareal_xbraid_errors.py . $fine_sim1 1

			mule.benchmark.cleanup_job_dirs || exit 1
		done

	fi


	echo ""
	echo ""
done


mule.benchmark.cleanup_all || exit 1

echo ""
echo_info "Test successful!"
