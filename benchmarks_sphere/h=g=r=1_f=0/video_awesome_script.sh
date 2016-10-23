#! /bin/bash

for i in linear_gaussian_*; do

	echo cp "$i/out_prog_h.mp4" "./video_$i.mp4"

done
