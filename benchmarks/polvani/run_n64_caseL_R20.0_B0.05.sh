#!/bin/bash
#SBATCH -o /home/martin/workspace/sweet/benchmarks/polvani/run_n64_caseL_R20.0_B0.05.txt
###SBATCH -e /home/martin/workspace/sweet/benchmarks/polvani/run_n64_caseL_R20.0_B0.05.txt.err
#SBATCH -D /home/martin/workspace/sweet/benchmarks/polvani
#SBATCH -J run_n64_caseL_R20.0_B0.05
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

mpiexec.hydra -genv OMP_NUM_THREADS 16 -ppn 1 ./build/swe_nonstaggered_advective_spectral_libfft_dealiasing_omp_intel_release  -X 1 -Y 1  -C 0.5 -R 4 -N 64 -g 399.999988079 -H 1 -S 1 -f 4.39822971503 -t 1000 -i "BINARY;../../tmp/polvani/CaseL/L64x640000.0000.h.raw;../../tmp/polvani/CaseL/L64x640000.0000.u.raw;../../tmp/polvani/CaseL/L64x640000.0000.v.raw" -v 2  -u 3.5527136788e-18 -U 8  -V 5  -O /home/martin/workspace/sweet/benchmarks/polvani/run_n64_caseL_R20.0_B0.05
