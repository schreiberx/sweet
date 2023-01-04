#! /bin/bash

cd "$(dirname $0)"

if [[ -z "$1" ]]; then
	TESTS=$(ls -1d ??_*/ ??_??_*/)
else
	TESTS=$@
fi

echo_info "Tests: ${TESTS}"

for TEST in $TESTS; do
	echo_info_hline
	echo_info "Executing script $i"
	echo_info_hline

	if [[ "${TEST: -7}" == "no_test" ]]; then
		echo_info "Skipping ${TEST} due to 'no_test' suffix"
		continue
	fi

	if [[ -e "$TEST/test.sh" ]]; then
		TEST_RUN="$TEST/test.sh"
	elif [[ -e "$TEST/test.py" ]]; then
		TEST_RUN="$TEST/test.py"
	else
		echo_error "test.sh or test.py not found for test '$TEST'"
		exit 1
	fi
	./$TEST_RUN || { echo "FAILED TO EXECUTE TEST $TEST_RUN"; exit 1; }

	echo_success_hline
	echo_success "Script $i successful"
	echo_success_hline
done

echo_success_hline
echo_success "*********************** FIN **************************"
echo_success_hline
