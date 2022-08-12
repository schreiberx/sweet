#! /bin/bash

for i in `find $@ -name "output_*.pdf"`; do
	if [[ "$i" == *"crop.pdf"* ]]; then
		continue
	fi

	echo "***********************************"
	echo "$i"
	echo "Cropping..."
	pdfcrop "$i"

	echo "Converting..."
	IN_FILE="${i/.pdf/-crop.pdf}"
	OUT_FILE="${i/.pdf/-crop.pdf}"
	echo "$IN_FILE -> $OUT_FILE"

	convert -density 600 "$IN_FILE" "$OUT_FILE"
done
