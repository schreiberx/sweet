#! /bin/bash

echo "**********************************************************"
echo "$0"
echo "**********************************************************"
echo "Arguments: $@"
echo "**********************************************************"

if [ -z "$2" ]; then
	echo "Usage: $0 [input files with wildcard] [output filename]"
	echo 
	echo "E.g.: $0 prog_h_ prog_h_output.mp4"
	exit 1
fi

if [ ! -z "$3" ]; then
	echo "Too many arguments"
	exit 1
fi

ffmpeg -y -framerate 30 -pattern_type glob -i "${1}" -c:v libx264 -pix_fmt yuv420p "${2}"
