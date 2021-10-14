#! /bin/bash

#if false; then
if true; then
	./pp_plot_csv.py output_prog_vrt_t00000000*csv || exit
	./pp_create_mp4.sh output_prog_vrt_ output_prog_vrt.mp4 || exit
fi


if false; then
#if true; then
	./pp_plot_csv.py output_prog_phi_t00000000*csv || exit
	./pp_create_mp4.sh output_prog_phi_ output_prog_phi.mp4 || exit
fi
