#! /bin/bash

if true; then

	# Which job directories to include in the data?
	DIRS=job_bench_*_C1.406e+01*
	DIRS+=" job_benchref_RT_tsm_ln_erk_tso4_tsob4_C2.000e+00_S864000"

	# Which tag in pickled job data to use to generate spectrum plots?
	SPECTAG="output_prog_h_pert_t00000864000.00000000_spectrum"

	# Output plot file
	OUTPUTFILE="plot_prog_h_pert_spectrum_simtime_10_days_C1.406e+01.pdf"

	# Plot!
	../pp_plot_spectrum.py $SPECTAG $OUTPUTFILE $DIRS

fi


if true; then

	# Which job directories to include in the data?
	DIRS=job_bench_*_C1.125e+02*
	DIRS+=" job_benchref_RT_tsm_ln_erk_tso4_tsob4_C2.000e+00_S864000"

	# Which tag in pickled job data to use to generate spectrum plots?
	SPECTAG="output_prog_h_pert_t00000864000.00000000_spectrum"

	# Output plot file
	OUTPUTFILE="plot_prog_h_pert_spectrum_simtime_10days_C1.125e+02.pdf"

	# Plot!
	../pp_plot_spectrum.py $SPECTAG $OUTPUTFILE $DIRS

fi
