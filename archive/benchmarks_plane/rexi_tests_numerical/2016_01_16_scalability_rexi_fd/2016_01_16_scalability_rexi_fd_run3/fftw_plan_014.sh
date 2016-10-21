#! /bin/bash

#SBATCH -o fftw_plan_014.txt
#SBATCH -J fftw_plan_014
#SBATCH --get-user-env
#SBATCH --clusters=mpp2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH --time=00:05:00

#declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS=1


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

cd /home/hpc/pr63so/di69fol/workspace/SWEET_2016_01_16/benchmarks_performance/rexi_tests/2016_01_06_scalability_rexi_fd_run3
cd ../../../

. local_software/env_vars.sh



mpiexec.hydra -genv OMP_NUM_THREADS 14 -envall -ppn 1 ./fftw_gen_wisdoms_all.sh 14 FFTW_WISDOM_nofreq_T14
