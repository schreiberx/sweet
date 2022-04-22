#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gnu
# Job id: RT_brossby_haurwitz_wave_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_l_irk_na_sl_settls_uv_only_tso2_tsob2_dt00005.00_W-00001_REXIDIRECT_M0032_N-001_X40031555.89280872_ptsm_l_irk_na_sl_settls_uv_only_pDt_-1_pStore_1
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/jcaldass/Development/IME/sweet"
source ./local_software/env_vars.sh "default_gnu" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


SCONS="scons  --mode=debug --compiler=gnu --sanitize= --debug-symbols=enable --simd=enable --mic=disable --fortran-source=disable --lapack=enable --program-binary-name= --sweet-mpi=disable --threading=omp --rexi-thread-parallel-sum=disable --benchmark-timings=disable --rexi-timings-additional-barriers=disable --rexi-allreduce=disable --program=swe_sphere --parareal=none --libpfasst=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=enable --plane-spectral-dealiasing=enable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --gui=disable --quadmath=disable -j 4"
echo "$SCONS"
$SCONS || exit 1


# chdir to execution directory
cd "/home/jcaldass/Development/IME/sweet/tests/70_program_swe_sphere_parareal/job_bench_COMP_plspec_pldeal_spspec_spdeal_fft_gnu_thomp_debug_RT_brossby_haurwitz_wave_g09.81_h010000.000_f1._8aea62dcf8be40df3b97085fe7275a67"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gnu
# Job id: RT_brossby_haurwitz_wave_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_l_irk_na_sl_settls_uv_only_tso2_tsob2_dt00005.00_W-00001_REXIDIRECT_M0032_N-001_X40031555.89280872_ptsm_l_irk_na_sl_settls_uv_only_pDt_-1_pStore_1
#


# mpiexec ... would be here without a line break
EXEC="/home/jcaldass/Development/IME/sweet/build/swe_sphere_COMP_plspec_pldeal_spspec_spdeal_fft_gnu_thomp_debug  -G 0 -g 9.80616 -H 10000 -f 0.00014584 -M 32 -N -1 --space-grid-use-c-staggering=0 -S 1 -X 40031555.89280872 -Y 40031555.89280872 --benchmark-name=rossby_haurwitz_wave -v 3 --dt=5.0 -o 25.0 -u 0.0 -t 500.0 --max-wallclock-time -1 -d 12 --timestepping-method=l_irk_na_sl_settls_uv_only --timestepping-order=2 --timestepping-order2=2 --rexi-method=direct --exp-direct-precompute-phin=0 --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=0 --reuse-plans=quick"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: default_gnu
# Job id: RT_brossby_haurwitz_wave_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_l_irk_na_sl_settls_uv_only_tso2_tsob2_dt00005.00_W-00001_REXIDIRECT_M0032_N-001_X40031555.89280872_ptsm_l_irk_na_sl_settls_uv_only_pDt_-1_pStore_1
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0