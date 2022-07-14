#! /bin/bash

cd "$(dirname $0)"

if [[ -z "$1" ]]; then
	TESTS=$(ls -1 ??_*/test.sh ??_*/test.py)
else
	TESTS=$@
fi

echo_info "Tests: ${TESTS}"

for i in $TESTS; do
	echo_info_hline
	echo_info "Executing script $i"
	echo_info_hline

	./$i || { echo "FAILED TO EXECUTE TEST $i"; exit 1; }

	echo_success_hline
	echo_success "Script $i successful"
	echo_success_hline
done

echo_success_hline
echo_success "*********************** FIN **************************"
echo_success_hline
