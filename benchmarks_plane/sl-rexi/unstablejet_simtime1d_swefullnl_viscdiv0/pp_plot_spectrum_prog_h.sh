#! /bin/bash

if true; then

	# Which job directories to include in the data?
	DIRS=job_bench_*_C1.406e+01*
	DIRS+=" job_benchref_RT_tsm_ln_erk_tso4_tsob4_C2.000e+00_S086400"

	# Which tag in pickled job data to use to generate spectrum plots?
	SPECTAG="output_prog_h_pert_t00000086400.00000000_spectrum"

	# Output plot file
	OUTPUTFILE="plot_prog_h_pert_simtime_1day.pdf"

	# Plot!
	../pp_plot_spectrum.py $SPECTAG $OUTPUTFILE $DIRS

fi
