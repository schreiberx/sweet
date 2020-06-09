#! /bin/bash

cd "$(dirname $0)"

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
	echo_info "Cleaning up test $i"
	echo_info_hline

	cd "$i" || exit

	rm -rf job_* || { echo "FAILED: $i"; exit 1; }
	rm -f output_* || { echo "FAILED: $i"; exit 1; }
	rm -f ./compile_platform_default_gnu.sh || { echo "FAILED: $i"; exit 1; }

	echo_success_hline
	echo_success "Cleanup of test $i successful"
	echo_success_hline
done

echo_success_hline
echo_success "*********************** FIN **************************"
echo_success_hline
