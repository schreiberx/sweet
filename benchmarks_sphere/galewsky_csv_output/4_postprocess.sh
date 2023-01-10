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

	if ! type ffmpeg >/dev/null 2>&1; then
		echo "''ffmpeg' required, but not found. Skipping video generation."
		exit 1
	else
		E=./pp_create_mp4.sh
		P1="${DIR}/output_prog_${FIELD}_t*.png"
	       	P2="${DIR}_output_prog_${FIELD}.mp4"
		A=($E "$P1" "$P2")
	
		echo "*****************************************************************************"
		echo "${A[@]}"
		echo "*****************************************************************************"

		"${A[@]}" || exit 1
	fi
done

