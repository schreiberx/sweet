#!/bin/bash

###############
## Unit tests following Section 3.6.2. "Debugging XBraid" from XBraid developer manual
## --> max_levels = 1 -> solution should be equal to the serial one
## --> max_levels = 1 and 2 processors in time -> idem
###############


get_tsm(){
	tsm=$1
	tsm=${tsm#*TS_};
	tsm=${tsm%.hpp};
	echo "$tsm";
}

dirname_serial="output_serial";
dirname_offline_error="output_simulations_offline_error";

cd "$(dirname $0)"

mule.benchmark.cleanup_job_dirs || exit 1
if [ -d $dirname_serial ]; then
	rm -rf $dirname_serial;
fi
mkdir $dirname_serial;

if [ -d $dirname_offline_error ]; then
	rm -rf $dirname_offline_error;
fi
mkdir $dirname_offline_error;



tsm_fine="dummy";
tsm_coarse="dummy";


## create and run reference simulation
echo "Running serial simulation"
./benchmarks_create.py ref $tsm_fine $tsm_coarse 1 > tmp_job_benchmark_create_dummy.txt || exit 1
mule.benchmark.jobs_run_directly || exit 1
mv job_bench_* "$dirname_serial"/.

# Get job directory name for this reference solution
# This will be reused throughout all other test cases
fine_sim=$(cat tmp_fine_sim.txt);

mule.benchmark.cleanup_job_dirs || exit 1

echo ""

echo "Running single level MGRIT simulations..."
for nproc in {1,2}; do

	echo "  -- nproc:" $nproc

	./benchmarks_create.py xbraid $tsm_fine $tsm_coarse $nproc > tmp_job_benchmark_create_dummy.txt || exit 1
	mule.benchmark.jobs_run_directly || exit 1
	cp -r "$dirname_serial"/"$fine_sim" .
	./compare_to_fine_solution.py $fine_sim;
	mule.benchmark.cleanup_job_dirs || exit 1
	echo ""
done;
