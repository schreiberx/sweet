#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gcc
# Job id: RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_dt5.00000000e-03_W-00001_REXIDIRECT_M0032_X40031555.89280872_plansquick
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/martin/workspace/sweet"
source ./local_software/env_vars.sh "default_gcc" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


SCONS="scons  --mode=debug --debug-symbols=enable --simd=enable --fortran-source=disable --lapack=disable --program-binary-name= --sweet-mpi=enable --threading=off --rexi-thread-parallel-sum=disable --benchmark-timings=disable --rexi-timings-additional-barriers=disable --rexi-allreduce=disable --program=parareal_ode --parareal=none --parareal-scalar=disable --parareal-plane=disable --parareal-sphere=disable --parareal-plane-swe=disable --parareal-plane-burgers=disable --xbraid=none --xbraid-scalar=disable --xbraid-plane=disable --xbraid-sphere=disable --xbraid-plane-swe=disable --xbraid-plane-burgers=disable --libpfasst=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=enable --plane-spectral-dealiasing=enable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --gui=disable --quadmath=disable -j 4"
echo "$SCONS"
$SCONS || exit 1


# chdir to execution directory
cd "/home/martin/workspace/sweet/tests/80_10_program_ode_scalar_GATHERED/80_20_program_ode_scalar_xbraid/job_bench_COMP_plspec_pldeal_spspec_spdeal_fft_mpi_debug_RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_55bdaedafaaebc83a7a7544e46289b26"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gcc
# Job id: RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_dt5.00000000e-03_W-00001_REXIDIRECT_M0032_X40031555.89280872_plansquick
#


# mpiexec ... would be here without a line break
EXEC="mpirun -n 1 /home/martin/workspace/sweet/build/parareal_ode_COMP_plspec_pldeal_spspec_spdeal_fft_mpi_debug  -G 0 -g 9.80616 -H 10000 -f 0.00014584 -M 32 --space-grid-use-c-staggering=0 -S 1 -X 40031555.89280872 -Y 40031555.89280872 --benchmark-name=unstablejet -v 3 --dt=0.005 -o 0.1 -u 0.0 -t 1.0 --max-wallclock-time -1 -d 12 --rexi-method=direct --exp-direct-precompute-phin=0 --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=1 --reuse-plans=quick"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: default_gcc
# Job id: RT_bunstablejet_g09.81_h010000.000_f1.458400e-04_u0.0_dt5.00000000e-03_W-00001_REXIDIRECT_M0032_X40031555.89280872_plansquick
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0