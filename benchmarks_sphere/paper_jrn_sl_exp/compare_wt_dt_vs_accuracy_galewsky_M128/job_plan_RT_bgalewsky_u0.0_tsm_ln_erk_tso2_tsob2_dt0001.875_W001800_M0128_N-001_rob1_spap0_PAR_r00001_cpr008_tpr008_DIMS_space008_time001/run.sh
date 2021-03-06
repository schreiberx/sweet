#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gnu
# Job id: RT_bgalewsky_u0.0_tsm_ln_erk_tso2_tsob2_dt0001.875_W001800_M0128_N-001_rob1_spap0
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/martin/workspace/sweet"
source ./local_software/env_vars.sh "default_gnu" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


# chdir to execution directory
cd "/home/martin/workspace/sweet/benchmarks_sphere/paper_jrn_sl_exp/compare_wt_dt_vs_accuracy_galewsky_M128/job_plan_RT_bgalewsky_u0.0_tsm_ln_erk_tso2_tsob2_dt0001.875_W001800_M0128_N-001_rob1_spap0_PAR_r00001_cpr008_tpr008_DIMS_space008_time001"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gnu
# Job id: RT_bgalewsky_u0.0_tsm_ln_erk_tso2_tsob2_dt0001.875_W001800_M0128_N-001_rob1_spap0
#


# mpiexec ... would be here without a line break
EXEC="/home/martin/workspace/sweet/build/swe_sphere_COMP_spspec_spdeal_numa2_fft_gnu_benchtime_mpi_thomp_release  -M 128 -N -1 --space-grid-use-c-staggering=0 -S 1 --benchmark-name=galewsky -v 0 --dt=1.875 -o -1 --output-file-name=- -u 0.0 -t 432000 --max-wallclock-time 1800 --instability-checks=0 -d 12 --timestepping-method=ln_erk --timestepping-order=2 --timestepping-order2=2 --rexi-method= --use-robert-functions=1 --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=0 --reuse-plans=-1"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: default_gnu
# Job id: RT_bgalewsky_u0.0_tsm_ln_erk_tso2_tsob2_dt0001.875_W001800_M0128_N-001_rob1_spap0
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0