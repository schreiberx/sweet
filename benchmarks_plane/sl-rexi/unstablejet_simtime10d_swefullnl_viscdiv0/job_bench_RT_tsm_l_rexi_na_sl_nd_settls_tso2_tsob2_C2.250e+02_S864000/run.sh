#! /bin/bash
#SBATCH -o /home/hpc/pr63qi/di69fol/workspace/sweet/benchmarks_plane/sl-rexi/unstablejet_simtime10d_swefullnl_viscdiv0/job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C2.250e+02_S864000/output.out
#SBATCH -e /home/hpc/pr63qi/di69fol/workspace/sweet/benchmarks_plane/sl-rexi/unstablejet_simtime10d_swefullnl_viscdiv0/job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C2.250e+02_S864000/output.err
#SBATCH -D /home/hpc/pr63qi/di69fol/workspace/sweet/benchmarks_plane/sl-rexi/unstablejet_simtime10d_swefullnl_viscdiv0/job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C2.250e+02_S864000
#SBATCH -J RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C2.250e+02_S864000
#SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
# the above is a good match for the
# CooLMUC2 architecture.
#SBATCH --mail-type=end 
#SBATCH --mail-user=schreiberx@gmail.com
#SBATCH --export=NONE 
#SBATCH --time=48:00:00

source /etc/profile.d/modules.sh


export OMP_NUM_THREADS=28
# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/hpc/pr63qi/di69fol/workspace/sweet"
source ./local_software/env_vars.sh "coolmuc_mpp2_gnu" || exit 1

export MPICC=mpigcc
export MPICXX=mpigxx
export MPIF90=mpifc


# chdir to execution directory
cd "/home/hpc/pr63qi/di69fol/workspace/sweet/benchmarks_plane/sl-rexi/unstablejet_simtime10d_swefullnl_viscdiv0/job_bench_RT_tsm_l_rexi_na_sl_nd_settls_tso2_tsob2_C2.250e+02_S864000"

# %SCRIPT_EXEC_PREFIX%


# mpiexec ... would be here without a line break
EXEC="/home/hpc/pr63qi/di69fol/workspace/sweet/build/swe_plane_COMP_plspec_pldeal_numa2_fft_gnu_thomp_release"
PARAMS=" -G 0 -g 9.80616 -H 10000 -f 0.00014584 -M 512 -N -1 --space-grid-use-c-staggering=0 -S 1 -X 40031555.89280872 -Y 40031555.89280872 --benchmark-name=unstablejet -v 3 --dt=225.0 -o 864000 -u 0.0 -t 864000 --instability-checks=0 -d 12 --timestepping-method=l_rexi_na_sl_nd_settls --timestepping-order=2 --timestepping-order2=2 --rexi-method=direct --use-robert-functions=1 --compute-error=1 --reuse-plans=-1"
echo "${EXEC} ${PARAMS}"

mpiexec -n 1 --perhost 1 $EXEC $PARAMS || exit 1

# %SCRIPT_EXEC_SUFFIX%
# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0