#! /bin/bash

#SBATCH -o /home/hpc/pr63so/di69fol/workspace/SWEET_2015_12_26/benchmarks_performance/rexi_tests_lrz_freq_waves/2015_12_27_scalability_rexi_spec/run_nr_fd_cart_cgrid_t008_n0032_r0001_a1.txt
###SBATCH -e /home/hpc/pr63so/di69fol/workspace/SWEET_2015_12_26/benchmarks_performance/rexi_tests_lrz_freq_waves/2015_12_27_scalability_rexi_spec/run_nr_fd_cart_cgrid_t008_n0032_r0001_a1.err
#SBATCH -J nr_fd_cart_cgrid_t008_n0032_r0001_a1
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH --time=01:00:00

#declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS=8


echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo


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

cd /home/hpc/pr63so/di69fol/workspace/SWEET_2015_12_26/benchmarks_performance/rexi_tests_lrz_freq_waves/2015_12_27_scalability_rexi_spec
cd ../../../

. local_software/env_vars.sh



# force to use FFTW WISDOM data
declare -x SWEET_FFTW_LOAD_WISDOM_FROM_FILE="FFTW_WISDOM_nofreq_T8"

time -p mpiexec.hydra -genv OMP_NUM_THREADS 8 -envall -ppn 3 -n 1 ./build/nr_fd_cart_cgrid_tyes_a1  --initial-freq-x-mul=2.0 --initial-freq-y-mul=1.0 -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -t 50 -R 4 -C 0.3 -N 32 -U 0 -S 0 --timestepping-mode 0 --staggering 1 -C 0.3

