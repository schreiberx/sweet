#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gnu
# Job id: RT_bgalewsky_linearbalance_u0.0_tsm_l_erk_tso2_tsob2_dt00001.88_W-00001_REXIDIRECT_M0128_N-001_plansquick
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/martin/workspace/sweet"
source ./local_software/env_vars.sh "default_gnu" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


# chdir to execution directory
cd "/home/martin/repositories/sweet/benchmarks_sphere/study_l_exp_special/compare_wt_dt_vs_accuracy_l_exp_special/job_plan_RT_bgalewsky_linearbalance_u0.0_tsm_l_erk_tso2_tsob2_dt00001.88_W-00001_REXIDIRECT_M0128_N-001_plansq_4a75e3bc225e3e05f714b5e50b213521"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gnu
# Job id: RT_bgalewsky_linearbalance_u0.0_tsm_l_erk_tso2_tsob2_dt00001.88_W-00001_REXIDIRECT_M0128_N-001_plansquick
#


# mpiexec ... would be here without a line break
EXEC="/home/martin/workspace/sweet/build/swe_sphere_COMP_spspec_spdeal_quadmath_fft_gnu_benchtime_thomp_release  -M 128 -N -1 --space-grid-use-c-staggering=0 -S 1 --benchmark-name=galewsky_linearbalance -v 0 --dt=1.875 -o -1 --output-file-name=- --output-file-mode=bin -u 0.0 -t 3600 --max-wallclock-time -1 --instability-checks=0 -d 12 --timestepping-method=l_erk --timestepping-order=2 --timestepping-order2=2 --rexi-method=direct --exp-direct-precompute-phin=0 --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=0 --reuse-plans=quick"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: default_gnu
# Job id: RT_bgalewsky_linearbalance_u0.0_tsm_l_erk_tso2_tsob2_dt00001.88_W-00001_REXIDIRECT_M0128_N-001_plansquick
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0