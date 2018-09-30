#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: martinium
# Job id: g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_ln_erk_tso2_tsob2_C000010_M0128
#




# Provide platform ID helper for scripts.
# This makes platform detection on the compute nodes easier since they might have mainly numerical identifiers
export SWEET_PLATFORM="martinium"

# Loading SWEET environment variables
cd "/home/martin/workspace/sweet"
source ./local_software/env_vars.sh "/home/martin/workspace/sweet/platforms/80_martinium/env_vars.sh" || exit 1



# chdir to execution directory
cd "/home/martin/workspace/sweet/benchmarks_sphere/paper_sph-rexi-nonlinear_compare_T_and_WT_vs_accuracy_galewsky_coolmuc/script_ln2g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_ln_erk_tso2_tsob2_C000010_M0128_MPI_space002_time001"



#
# Generating function: jobscript_get_exec_prefix
# Platform: martinium
# Job id: g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_ln_erk_tso2_tsob2_C000010_M0128
#




# SHTNS plans required which have to be precomputed once.
# You can generate them by simply running one benchmark.
# Then, copy them to the main benchmark folder
if [ -e ../shtns_cfg ]; then
	ln ../shtns_cfg ./ -sf || exit 1
	ln ../shtns_cfg_fftw ./ -sf || exit 1
fi

echo "CPU Frequencies:"
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq | sort -u
echo ""

# Manual override in jobs_create.py
# Override OMP_NUM_THREADS and use MASTER binding
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=MASTER



#
# Generating function: jobscript_get_exec_command
# Platform: martinium
# Job id: g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_ln_erk_tso2_tsob2_C000010_M0128
#


# mpiexec ... would be here without a line break
EXEC="/home/martin/workspace/sweet/build/swe_sphere_spspec_spdeal_quadmath_mpi_rxtime_numa2_fft_gnu_release  -g 9.81 -H 10000 -f 7.2921e-05 -F 0 -a 6371220 -M 128 --pde-id 0 --staggering=0 -S 1 -X 1 -Y 1 -s 4 --benchmark= -v 0 --dt=10 -o 432000 -u 0.0 -t 432000 --stability-checks=1 --use-linear-div=0 --use-local-visc=0 --timestepping-method=ln_erk --timestepping-order=2 --timestepping-order2=2 --normal-mode-analysis-generation=0 --rexi-method= --polvani-rossby=-1.0 --polvani-froude=-1.0 --use-robert-functions=1 --compute-error=0 --shtns-use-plans=1"
echo "$EXEC"
$EXEC



echo "CPU Frequencies:"
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq | sort -u
echo ""



#
# Generating function: jobscript_get_exec_suffix
# Platform: martinium
# Job id: g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_ln_erk_tso2_tsob2_C000010_M0128
#




#
# Generating function: jobscript_get_footer
# Platform: martinium
# Job id: g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_ln_erk_tso2_tsob2_C000010_M0128
#


