#! /usr/bin/env bash

#
# Clean junk from output files, keep only MULE lines
#


if [[ -z "$1" ]]; then
	DIRS=bench*
else
	DIRS=$@
fi

P="$(pwd)"
echo $P
for JOBDIR in $DIRS; do
	cd "$JOBDIR" || exit 1
    echo $JOBDIR
	SUBDIRS=job_bench*
	
	for DIR in $SUBDIRS; do
		echo $DIR
		out=$DIR/output.out
		
		out_orig=$DIR/output_orig.out
		out_clean=$DIR/output_clean.out

		#Safety checks
		size="$(wc -l <"$out")"
		if test -f "$out_orig"; then
			size_orig="$(wc -l <"$out_orig")"
			if [ $size_orig > $size ]; then
				echo "   Cleanning already happened here"
			fi
		else
			cp $out $out_clean
			cp $out $out_orig
		fi

		#Clear out file
		#head $out_orig
		grep "MULE" $out_orig > $out_clean 
		cp $out_clean $out

		#Safety checks
		size_orig="$(wc -l <"$out_orig")"
		size_clean="$(wc -l <"$out_clean")"
		echo "   Sizes (out, original, clean): "  $size, $size_orig, $size_clean
	done
	cd "$P"
done