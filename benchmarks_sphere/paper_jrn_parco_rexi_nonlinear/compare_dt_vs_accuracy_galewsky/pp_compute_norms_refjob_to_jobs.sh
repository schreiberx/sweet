#! /usr/bin/env bash

#
# Compute error norms between
#
# one reference job (directory as 1st argument)
#
# other jobs (2nd to last argument)
#
# The output data is determined by the reference jobs files
#	output_*.csv
#
# A filter "0\{11\}.0\{8\}" is applied to sort out the 0-time stamp values
#



if [[ -z "$3" ]]; then
	echo ""
	echo "	Usage:"
	echo "		$0 [reference job directory] [reference data ending] [job directory 1] [job directory 2] ..."
	echo ""
	echo "	Example:"
	echo "		$0 ./*_ref_* 00000000120.00000000.csv ./job_jadda_jadda1 ./job_foo_bar2"
	echo ""
fi

REF_JOB="$1"
REF_FILE_ENDING="$2"

CMP_JOBS="${@:3}"

DIFF_FILENAMES="$( cd $REF_JOB; ls -1 output_*.csv | grep "${REF_FILE_ENDING}\$")"
echo_info_hline
echo_info "Using reference files:"
echo_info "$DIFF_FILENAMES"
echo_info_hline

for DIFF_FILENAME in $DIFF_FILENAMES; do
	echo "Comparison filename: ${DIFF_FILENAME}"

	REF_FILE="${REF_JOB}/${DIFF_FILENAME}"
	REF_TAG="$(echo $DIFF_FILENAME | sed "s/output_//" | sed "s/_t.*//")"

	if [[ -z "$REF_TAG" ]]; then
		echo_error_exit "Invalid reference tag!"
	fi


	for CMP_JOB in $CMP_JOBS; do
		if [[ "$CMP_JOB" = "$REF_JOB" ]]; then
			continue
		fi

		CMP_FILE="${CMP_JOB}/${DIFF_FILENAME}"
		CMP_PICKLE="${CMP_JOB}/sphere_data_diff_${REF_TAG}.pickle"

		echo ""
		echo "Computing differences between"
		echo " FileA: ${REF_FILE}"
		echo " FileB: ${CMP_FILE}"
#		echo " Reference tagname: ${REF_TAG}"
		echo " Out pickle file: ${CMP_PICKLE}"
		echo ""

		if [[ ! -e "$CMP_FILE" ]]; then
			echo "Skipping this job, since output doesn't exist (very likely due to an instability)"
			continue
		fi

		# Ignore errors since they are likely to indicate non-existing output data due to instabilities
		echo_exec ../../../python_mods/SphereDataDiff.py "$REF_FILE" "$CMP_FILE" "$CMP_PICKLE"
	done

done

echo "FIN"
