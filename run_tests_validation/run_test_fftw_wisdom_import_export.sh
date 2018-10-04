#! /bin/bash

echo_info_hline
echo_info "Running tests for REXI"
echo_info_hline

cd ..

scons --unit-test=test_fftw_wisdom_import_export --gui=disable --mode=release || echo_error_exit "Failed scons"

RES=128

echo_info_hline
echo_info_hline
echo_info_hline
echo_info "This unit test checks if the Wisdom load/store works properly"
echo_info "There's an issue on 'martinium' (2 cores, 4 threads) if using"
echo_info "more than 2 threads for the FFTW plans. This results in no"
echo_info "plans being generated"
echo_info_hline
echo_info_hline
echo_info_hline

# Determine number of processors
NPROC=$(nproc --all)

# Use twice as many processors
NPROC=$(($NPROC*2))

# Generate testing sequence 
SEQ=$(seq ${NPROC})

# Make pretty
SEQ=$(echo -n "${SEQ}" | tr '\n' ' ')

# Add autodetection
SEQ="-1 ${SEQ}"

# Output
echo_info " Testing for THREADS \in { ${SEQ} }"

for THREADS in $SEQ; do

	echo_info "Cleaning up file 'sweet_fftw'"
	rm -rf "sweet_fftw"

	echo_info_hline
	echo_info "Testing for '${THREADS}' threads"
	echo_info_hline

	# Generate wisdom
	EXEC="./build/test_fftw_*_release 128 128 $THREADS 1"
	echo "$EXEC"
	$EXEC || echo_error_exit "FAILED generation of wisdom for '$THREADS' threads"

	# Reuse wisdom and stop if 
	EXEC="./build/test_fftw_*_release 128 128 $THREADS 2"
	echo "$EXEC"
	$EXEC || echo_error_exit "FAILED reutilization of wisdom for '$THREADS' threads"


	echo_success_hline
	echo_success "Wisdom test successful for '$THREADS' threads"
	echo_success_hline

done


echo_success_hline
echo_success " All tests for THREADS \in { $SEQ } passed"
echo_success_hline
