#! /bin/bash

if [ -z "$2" ]; then
	echo "Usage: $0 [prefix_of_filename] [output filename]"
	echo 
	echo "E.g.: $0 prog_h_ prog_h_output.mp4"
	exit
fi

ffmpeg -y -framerate 30 -pattern_type glob -i "$1*.png" -c:v libx264 -pix_fmt yuv420p $2
