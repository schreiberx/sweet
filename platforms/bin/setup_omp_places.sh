#! /usr/bin/env bash

if [[ -z "$2" ]]; then
	echo "Usage:"
	echo "	./$0 [oversubscription|nooversubscription] [close|spread] ..."
	return
fi

function join()
{
	local IFS="$1";
	shift;
	echo "$*";
}

DIR=$(dirname $0)

if [[ "$1" = "nooversubscription" ]]; then

	if [ "$2" = "close" ]; then
		CORE_IDS=$(${SWEET_ROOT}/platforms/bin/get_physical_cpus_close.sh)
	elif [ "$2" = "spread" ]; then
		echo "TODO: Not yet supported"
		return
	else
		echo "Unknown placement '$2'"
		return
	fi

elif [[ "$1" = "oversubscription" ]]; then
	echo "TODO: oversubscription"
	return
fi

NUM_CORES="$(echo "$CORE_IDS" | wc -l)"

LIST="$(join , $CORE_IDS)"
LIST="{$(echo $LIST | sed "s/,/},{/g")}"

if [[ -z "$LIST" ]]; then
	echo "$CORE_IDS"
	echo "Unknown error"
	return
fi

#export OMP_NUM_THREADS="$NUM_CORES"
#echo "export OMP_NUM_THREADS=$NUM_CORES"

export OMP_PLACES="$LIST"
echo "export OMP_PLACES=$OMP_PLACES"
