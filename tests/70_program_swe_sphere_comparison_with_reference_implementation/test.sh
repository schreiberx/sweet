#! /bin/bash


echo_info_hline
echo_info "Running comparisons of SWEET SWE on the sphere implementation with Python reference implementation"
echo_info_hline

cd "$(dirname $0)"


echo_info "Cleaning up..."
mule.benchmark.cleanup_all || exit 1

echo_info "Generating benchmark job scripts..."
./benchmark_create_job_scripts || exit 1

echo_info "Generating reference solution..."
./benchmark_gen_swe_reference_solution.py || exit 1

echo_info "Running various time steppers..."
mule.benchmark.jobs_run_directly || exit 1

echo_info "Postprocessing..."
./postprocessing || exit 1

echo_info "Cleaning up..."
mule.benchmark.cleanup_all || exit 1

echo_success "Success!"
