#! /bin/bash


. /etc/profile.d/modules.sh

module unload gcc
module unload fftw

module unload python
module load python/2.7_anaconda_nompi


module unload intel
module load intel/16.0

module unload mpi.intel
module load mpi.intel/5.1

module load gcc/5

cd /home/martin/workspace/sweet/benchmarks/rexi_tests_lrz/2016_10_16_scalability_nr_fd
cd ../../../

. local_software/env_vars.sh



scons --compiler=intel --sweet-mpi=enable --program=swe_rexi --plane-spectral-space=disable --libfft=enable --plane-spectral-dealiasing=disable --mode=release --numa-block-allocator=1  --threading=omp --program-binary-name=nr_fd_cart_cgrid_tyes_a1 || exit 1
