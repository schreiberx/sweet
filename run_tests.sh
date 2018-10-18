#! /bin/bash

source ./local_software/env_vars.sh

# only execute scripts which do not start with an underscore
for i in tests/??_*/test.sh; do
	echo_info_hline
	echo_info "Executing script $i"
	echo_info_hline

	./$i || exit 1
done

echo_info_hline
echo_info "*********************** FIN **************************"
echo_info_hline
