#! /bin/bash

TAGS=()
TAGS+=("simulation_successfully_finished")
TAGS+=("DUE TO TIME LIMIT")
TAGS+=("INSTABILITY DETECTED")

#for TAG in "${TAGS[@]}"; do
#	echo "$TAG"
#done

# Really! We need to use "..." marks to hand over arrays

mule.benchmark.jobs_get_missing_tag "${TAGS[@]}"

