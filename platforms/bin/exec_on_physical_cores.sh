#! /usr/bin/env bash

if [[ -z "$1" ]]; then
	echo_info "Usage:"
	echo_info "	$0 [exec] [parameter1] [parameter2] ..."
	exit 1
fi

function join()
{
	local IFS="$1";
	shift;
	echo "$*";
}

DIR=$(dirname $0)
CORE_IDS=$(${DIR}/get_physical_cpus_close.sh)

LIST=$(join , $CORE_IDS)

if [[ -z "$LIST" ]]; then
	echo "$CORE_IDS"
	echo_error_exit "Unknown error"
fi

EXEC="taskset -c $LIST $@"
echo_info "Executing '${EXEC}'"
$EXEC || exit 1


