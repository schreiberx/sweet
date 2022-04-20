#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gnu
# Job id: RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_l_erk_tso2_tsob2_dt00005.00_W-00001_REXIDIRECT_M0008_N-001_X40031555.89280872_plansquick_par_3_ptsm_l_erk_pDt_-1_pStore_0
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/jcaldass/Development/IME/sweet"
source ./local_software/env_vars.sh "default_gnu" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


SCONS="scons  --mode=release --compiler=gnu --sanitize= --debug-symbols=enable --simd=enable --mic=disable --fortran-source=disable --lapack=enable --program-binary-name= --sweet-mpi=disable --threading=omp --rexi-thread-parallel-sum=disable --benchmark-timings=disable --rexi-timings-additional-barriers=disable --rexi-allreduce=disable --program=swe_plane --parareal=serial --libpfasst=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=enable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --gui=disable --quadmath=disable -j 4"
echo "$SCONS"
$SCONS || exit 1


# chdir to execution directory
cd "/home/jcaldass/Development/IME/sweet/tests/70_program_swe_plane_parareal/job_bench_COMP_plspec_spspec_spdeal_fft_gnu_thomp_release_RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0._9cd0ef0a7e252ab883612eb0e6ff0cfc"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gnu
# Job id: RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_l_erk_tso2_tsob2_dt00005.00_W-00001_REXIDIRECT_M0008_N-001_X40031555.89280872_plansquick_par_3_ptsm_l_erk_pDt_-1_pStore_0
#


# mpiexec ... would be here without a line break
EXEC="/home/jcaldass/Development/IME/sweet/build/swe_plane_COMP_plspec_spspec_spdeal_fft_gnu_thomp_release  -G 0 -g 9.80616 -H 10000 -f 0.00014584 -M 8 -N -1 --space-grid-use-c-staggering=0 -S 1 -X 40031555.89280872 -Y 40031555.89280872 --benchmark-name=unstablejet -v 3 --dt=5.0 -o 10.0 -u 0.0 -t 30.0 --max-wallclock-time -1 -d 12 --timestepping-method=l_erk --timestepping-order=2 --timestepping-order2=2 --rexi-method=direct --exp-direct-precompute-phin=0 --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=1 --reuse-plans=quick --parareal-enable=1 --parareal-coarse-slices=3 --parareal-convergence-threshold=-1 --parareal-verbosity=6 --parareal-max-simulation-time=30.0 --parareal-coarse-timestepping-method=l_erk --parareal-coarse-timestepping-order=2 --parareal-coarse-timestepping-order2=2 --parareal-coarse-timestep-size=-1 --parareal-load-ref-csv-files=1 --parareal-path-ref-csv-files=../simulations_offline_error/l_erk_l_erk/job_benchref_COMP_plspec_spspec_spdeal_fft_gnu_thomp_release_RT_bunstablejet_g09.81_h010000.000_f1.458400e-04__9f7b1932886c41796bf25ac88b630611 --parareal-load-fine-csv-files=1 --parareal-path-fine-csv-files=../simulations_offline_error/l_erk_l_erk/job_bench_COMP_plspec_spspec_spdeal_fft_gnu_thomp_release_RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0._ea3e7d284afbbea69a0d153f32dcda8e --parareal-store-iterations=0"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: default_gnu
# Job id: RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_tsm_l_erk_tso2_tsob2_dt00005.00_W-00001_REXIDIRECT_M0008_N-001_X40031555.89280872_plansquick_par_3_ptsm_l_erk_pDt_-1_pStore_0
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0