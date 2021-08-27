#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: ppeixoto_usp_gnu
# Job id: RT_bwilliamson2_linear_u0.0_tsm_ln_erk_tso2_tsob2_dt00007.50_W-00001_M0064_N-001_spap0_plansquick
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/media/pedrosp/Data/Simulations/sweet"
source ./local_software/env_vars.sh "ppeixoto_usp_gnu" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


# chdir to execution directory
cd "/media/pedrosp/Data/Simulations/sweet/benchmarks_sphere/ppeixoto/compare_wt_dt_vs_accuracy_galewsky_reprod_2020_05_30/job_plan_RT_tsm_ln_erk_dt00007.50_W-00001"


export OMP_NUM_THREADS=16
export OMP_DISPLAY_ENV=VERBOSE

echo "Affnity: compact"
source $MULE_ROOT/platforms/bin/setup_omp_places.sh nooversubscription close

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: ppeixoto_usp_gnu
# Job id: RT_bwilliamson2_linear_u0.0_tsm_ln_erk_tso2_tsob2_dt00007.50_W-00001_M0064_N-001_spap0_plansquick
#


# mpiexec ... would be here without a line break
EXEC="/media/pedrosp/Data/Simulations/sweet/build/swe_sphere_COMP_spspec_spdeal_quadmath_numa2_fft_gnu_benchtime_mpi_thomp_release"
PARAMS=" -M 64 -N -1 --space-grid-use-c-staggering=0 -S 1 --benchmark-name=williamson2_linear -v 0 --dt=7.5 -o -1 --output-file-name=- --output-file-mode=bin -u 0.0 -t 43200 --max-wallclock-time -1 --instability-checks=0 -d 12 --timestepping-method=ln_erk --timestepping-order=2 --timestepping-order2=2 --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=1 --reuse-plans=quick"
echo "${EXEC} ${PARAMS}"

mpiexec -n 1 $EXEC $PARAMS || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: ppeixoto_usp_gnu
# Job id: RT_bwilliamson2_linear_u0.0_tsm_ln_erk_tso2_tsob2_dt00007.50_W-00001_M0064_N-001_spap0_plansquick
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0