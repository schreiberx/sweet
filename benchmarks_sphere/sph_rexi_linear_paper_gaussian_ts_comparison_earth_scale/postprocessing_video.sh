#! /bin/bash

cd "output_video_rk4_ts050"
for i in output_prog_h_*csv; do
	../pp_plot_csv.py $i
done

../pp_create_mp4.sh output_prog_h_ ../output_video_rk4_ts050.mp4
