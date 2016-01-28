#!/bin/bash
#SBATCH -o /home/martin/workspace/sweet/benchmarks/polvani/run_n128_caseE_R0.25_B0.2.txt
###SBATCH -e /home/martin/workspace/sweet/benchmarks/polvani/run_n128_caseE_R0.25_B0.2.txt.err
#SBATCH -D /home/martin/workspace/sweet/benchmarks/polvani
#SBATCH -J run_n128_caseE_R0.25_B0.2
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=end
#SBATCH --mail-user=martin.schreiber@in.tum.de
#SBATCH --export=NONE
#SBATCH --time=24:00:00

source /etc/profile.d/modules.sh

module load fftw/mpi/3.3

. ~/bin/local_vars.sh
. ~/bin/intel_vars.sh


declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
cd ../../
rm -f FFTW_WISDOM.cached

scons --compiler=intel --program=swe_nonstaggered_advective --spectral-space=enable --libfft=enable --spectral-dealiasing=enable --threading=omp --mode=release

mpiexec.hydra -genv OMP_NUM_THREADS 16 -ppn 1 ./build/swe_nonstaggered_advective_spectral_libfft_dealiasing_omp_intel_release  -X 1 -Y 1  -C 0.5 -R 4 -N 128 -g 24.9999992549 -H 1 -S 1 -f 351.858377202 -t 1000 -i "BINARY;../../tmp/polvani/CaseE/E128x1280000.0000.h.raw;../../tmp/polvani/CaseE/E128x1280000.0000.u.raw;../../tmp/polvani/CaseE/E128x1280000.0000.v.raw" -v 2  -u 1.38777878078e-20 -U 8  -V 5  -O /home/martin/workspace/sweet/benchmarks/polvani/run_n128_caseE_R0.25_B0.2
