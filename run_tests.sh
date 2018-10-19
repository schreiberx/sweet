#! /bin/bash

source ./local_software/env_vars.sh

if [[ -z "$1" ]]; then
	TESTS=$(ls -1 tests/??_*/test.sh tests/??_*/test.py)
else
	TESTS=$@
fi

echo_info "Tests: ${TESTS}"

for i in $TESTS; do
	echo_info_hline
	echo_info "Executing script $i"
	echo_info_hline

	./$i || exit 1

	echo_success_hline
	echo_success "Script $i successful"
	echo_success_hline
done

echo_success_hline
echo_success "*********************** FIN **************************"
echo_success_hline
