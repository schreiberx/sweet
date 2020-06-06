#! /bin/bash

cd "$(dirname $0)"

export

if [[ -z "$1" ]]; then
	TESTS=$(ls -1d test_*)
else
	TESTS=$@
fi

echo_info "Running: ${TESTS}"
p="$(pwd)"

for i in $TESTS; do
	cd "$p"

	echo_info_hline
	echo_info "Executing script $i"
	echo_info_hline

	cd "$i" || exit
	./benchmark_create_jobs.py || { echo "FAILED: $i"; exit 1; }
	./compile_platform_default_gnu.sh || { echo "FAILED: $i"; exit 1; }
	mule.benchmark.jobs_run_directly || { echo "FAILED: $i"; exit 1; }
	./postprocessing.sh || { echo "FAILED: $i"; exit 1; }

	echo_success_hline
	echo_success "Script $i successful"
	echo_success_hline
done

echo_success_hline
echo_success "*********************** FIN **************************"
echo_success_hline
