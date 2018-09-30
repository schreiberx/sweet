#!/bin/bash
# TARGET MACHINE: cheyenne
#
## project code
#PBS -A NCIS0002
## regular limit: 12 hours
## economy queue
#PBS -q economy
## shared queue
######PBS -q share
## wall-clock time (hrs:mins:secs)
#PBS -l walltime=00:10:00
## select: number of nodes
## ncpus: number of CPUs per node
## mpiprocs: number of ranks per node
#PBS -l select=64:ncpus=36:mpiprocs=2:ompthreads=18
#PBS -l select=cpufreq=2300000
#
#PBS -N script_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_l_rexi_n_erk_ver1_tso2_tsob2_C000720_REXI
#PBS -o /gpfs/u/home/martins/workspace/sweet/benchmarks_sphere/paper_sph-rexi-nonlinear_compare_time_to_solution_vs_accuracy_galewsky_cheyenne/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_l_rexi_n_erk_ver1_tso2_tsob2_C000720_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128/output.out
#PBS -e /gpfs/u/home/martins/workspace/sweet/benchmarks_sphere/paper_sph-rexi-nonlinear_compare_time_to_solution_vs_accuracy_galewsky_cheyenne/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_l_rexi_n_erk_ver1_tso2_tsob2_C000720_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128/output.err

export OMP_NUM_THREADS=18





# Manual override in jobs_create.py
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=MASTER



cd "/gpfs/u/home/martins/workspace/sweet/benchmarks_sphere/paper_sph-rexi-nonlinear_compare_time_to_solution_vs_accuracy_galewsky_cheyenne/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_l_rexi_n_erk_ver1_tso2_tsob2_C000720_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128"

BASEDIR="`pwd`"

SWEETROOT="/gpfs/u/home/martins/workspace/sweet/benchmarks_sphere/paper_sph-rexi-nonlinear_compare_time_to_solution_vs_accuracy_galewsky_cheyenne/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_l_rexi_n_erk_ver1_tso2_tsob2_C000720_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time128/../../../"
cd "$SWEETROOT"
pwd

# Always load local software
#is this really the root?
if test -e ./local_software/env_vars.sh ; then
	source ./local_software/env_vars.sh || exit 1
else
	echo "Warning: changing SWEETROOT directory"	
	cd ..
	SWEETROOT="`pwd`"
	pwd
	source ./local_software/env_vars.sh || exit 1
fi

###make clean || exit 1


cd "$BASEDIR"
pwd



# SHTNS plans required which have to be precomputed once.
# You can generate them by simply running one benchmark.
# Then, copy them to the main benchmark folder
ln ../shtns_cfg ./ -sf || exit 1
ln ../shtns_cfg_fftw ./ -sf || exit 1

echo "CPU Frequencies:"
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq | sort -u
echo ""

EXEC="$SWEETROOT/build/swe_sphere_spspec_spdeal_quadmath_mpi_rxtime_numa2_fft_intel_release  -g 9.81 -H 10000 -f 7.2921e-05 -F 0 -a 6371220 -M 128 --pde-id 0 --staggering=0 -S 1 -X 1 -Y 1 -s 4 --benchmark= -v 0 --dt=720 -o 432000 -u 0.0 -t 432000 --stability-checks=1 --use-linear-div=0 --use-local-visc=0 --timestepping-method=l_rexi_n_erk_ver1 --timestepping-order=2 --timestepping-order2=2 --normal-mode-analysis-generation=0 --rexi-method=ci --rexi-half=0 --rexi-normalization=0 --rexi-sphere-preallocation=1 --rexi-use-direct-solution=0 --rexi-ext-modes=0 --rexi-ci-n=128 --rexi-ci-max-real=10.0 --rexi-ci-max-imag=30.0 --rexi-ci-primitive=circle --rexi-ci-gaussian-filter-scale=0.0 --rexi-ci-gaussian-filter-dt-norm=0.0 --rexi-ci-gaussian-filter-exp-N=0.0 --polvani-rossby=-1.0 --polvani-froude=-1.0 --use-robert-functions=1 --compute-error=0 --shtns-use-plans=0"


echo "$EXEC"
pwd
#ln -s "$SWEETROOT/data/" "$BASEDIR/data"   #Symlink for GUI directory, if necessary
mpiexec_mpt -n 128  omplace  -nt 18  -vv $EXEC || exit 1



echo "CPU Frequencies:"
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_cur_freq | sort -u
echo ""

