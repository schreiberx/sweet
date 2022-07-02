#! /bin/bash

FIELD=vrt
for DIR in job_*; do
	echo "${DIR}"

	if true; then
		# Faster version since SH transformations need to be set up only once
		./pp_plot.py ${DIR}/output_prog_${FIELD}_t*.csv || exit 1
	else
		for FILE in ${DIR}/output_prog_${FIELD}_t*.csv; do
			echo " + ${FILE}"
			./pp_plot.py "${FILE}" || exit 1
		done
	fi

	./pp_create_mp4.sh "${DIR}/output_prog_${FIELD}_t*.png" ${DIR}_output_prog_${FIELD}.mp4 || exit 1
done

