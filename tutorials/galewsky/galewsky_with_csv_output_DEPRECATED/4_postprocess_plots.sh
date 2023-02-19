#! /bin/bash

FIELD=vrt
EXT=csv
for DIR in job_*/; do
	echo "*****************************************************************************"
	echo "* ${DIR}"
	echo "*****************************************************************************"

	# Remove trailing slash
	DIR=${DIR%/}

	if true; then
		# Faster version since SH transformations need to be set up only once
		./pp_plot.py ${DIR}/output_prog_${FIELD}_t*.${EXT} || exit 1
	else
		for FILE in ${DIR}/output_prog_${FIELD}_t*.${EXT}; do
			echo " + ${FILE}"
			./pp_plot.py "${FILE}" || exit 1
		done
	fi
done

