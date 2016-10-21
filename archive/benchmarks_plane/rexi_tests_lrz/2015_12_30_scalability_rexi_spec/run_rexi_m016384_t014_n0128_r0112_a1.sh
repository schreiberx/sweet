#! /bin/bash

#SBATCH -o /home/hpc/pr63so/di69fol/workspace/SWEET_2015_12_26/benchmarks_performance/rexi_tests_lrz_freq_waves/2015_12_30_scalability_rexi_spec/run_rexi_m016384_t014_n0128_r0112_a1.txt
###SBATCH -e /home/hpc/pr63so/di69fol/workspace/SWEET_2015_12_26/benchmarks_performance/rexi_tests_lrz_freq_waves/2015_12_30_scalability_rexi_spec/run_rexi_m016384_t014_n0128_r0112_a1.err
#SBATCH -J rexi_m016384_t014_n0128_r0112_a1
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --ntasks=112
#SBATCH --cpus-per-task=14
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH --time=03:00:00

#declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS=14


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

cd /home/hpc/pr63so/di69fol/workspace/SWEET_2015_12_26/benchmarks_performance/rexi_tests_lrz_freq_waves/2015_12_30_scalability_rexi_spec
cd ../../../

. local_software/env_vars.sh



# force to use FFTW WISDOM data
declare -x SWEET_FFTW_LOAD_WISDOM_FROM_FILE="FFTW_WISDOM_nofreq_T14"

time -p mpiexec.hydra -genv OMP_NUM_THREADS 14 -envall -ppn 2 -n 112 ./build/rexi_m_tyes_a1  --initial-freq-x-mul=2.0 --initial-freq-y-mul=1.0 -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -t 50 -R 4 -C 0.3 -N 128 -U 0 -S 0 --use-specdiff-for-complex-array 1 --rexi-h 0.2 --timestepping-mode 1 --staggering 0 --rexi-m=16384 -C -5.0

