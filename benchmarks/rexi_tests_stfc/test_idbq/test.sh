#BSUB -o /gpfs/stfc/local/HCP011/smp03/mxs64-smp03/workspace/SWEET_2015_12_16/benchmarks/rexi_tests/2015_12_15_error_plot_fd/test.txt
#BSUB -e /gpfs/stfc/local/HCP011/smp03/mxs64-smp03/workspace/SWEET_2015_12_16/benchmarks/rexi_tests/2015_12_15_error_plot_fd/test.err
#BSUB -R "span[ptile=1]"
#BSUB -R "affinity[core(24):cpubind=core:distribute=pack]"
#BSUB -R same[type:model]
###BSUB -R order[hosts]
# All on the same rack - TODO: does this really work?
##BSUB -R "cu[type=rack]"
##BSUB -R "cu[maxcus=1]"
# dedicated network
###BSUB -network "type=sn_all: usage=dedicated"
###BSUB -network "type=sn_single: usage=dedicated"
# exclusive resource
#BSUB -x
##BSUB -m ib0
#BSUB -J rexi_fd_par_m0008_t024_n0032_r0001_a1
#BSUB -W 01:00
#BSUB -n 1

echo "LSB_BIND_CPU_LIST"
echo "$LSB_BIND_CPU_LIST"

echo "LSB_BIND_MEM_LIST"
echo "$LSB_BIND_MEM_LIST"

echo "LSB_BIND_MEM_POLICY"
echo "$LSB_BIND_MEM_POLICY"

#    RM_CPUTASKn"
echo "RM_MEM_AFFINITY"
echo "$RM_MEM_AFFINITY"

echo "OMP_NUM_THREADS"
echo "$OMP_NUM_THREADS"

echo



. /etc/profile.d/modules.sh

#module load gcc/5.1.0
module unload gcc

# Manual GCC override, since 'module load gcc' is not possible after loading intel module
#export INCLUDE_PATH="/gpfs/stfc/local/apps/gcc/gcc/4.9.1/include:$INCLUDE_PATH"
#export LD_LIBRARY_PATH="/gpfs/stfc/local/apps/gcc/gcc/4.9.1/lib64:/gpfs/stfc/local/apps/gcc/gcc/4.9.1/lib:/gpfs/stfc/local/apps/gcc/utilities/lib:$LD_LIBRARY_PATH"
#export PATH="/gpfs/stfc/local/apps/gcc/gcc/4.9.1/bin:/gpfs/stfc/local/apps/gcc/utilities/bin:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/opt/xcat/bin:/opt/xcat/sbin:/opt/xcat/share/xcat/tools:/usr/local/bin:/usr/local/X11/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/usr/lpp/mmfs/bin:$PATH"


# This FFTW does not support OMP!!!
#module load fftw/3.3.4
module load python/2.7.8

#module load gcc/5.1.0
#module unload intel/license
module load intel/15.2.164
module load intel_mpi

cd /gpfs/stfc/local/HCP011/smp03/mxs64-smp03/workspace/SWEET_2015_12_16/benchmarks/rexi_tests/2015_12_15_error_plot_fd
cd ../../../

. ~/bin/local_vars.sh
. ~/bin/intel_vars.sh

. local_software/env_vars.sh


#declare -x NUMA_BLOCK_ALLOC_VERBOSITY=1
# this ignored hyperthreads
declare -x KMP_AFFINITY="granularity=thread,compact,1,0"
declare -x OMP_NUM_THREADS=24

# force to use FFTW WISDOM data
declare -x SWEET_FFTW_LOAD_WISDOM_FROM_FILE="FFTW_WISDOM_T0"

./build/rexi_fd_par_m_tno_a1  -s 5 -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -t 50 -R 4 -C 0.3 -N 32 -U 0 -S 0 --use-fd-for-complex-array 1 --rexi-h 0.8 --timestepping-mode 1 --staggering 0 --rexi-m=8 -C -5.0 #> /gpfs/stfc/local/HCP011/smp03/mxs64-smp03/workspace/SWEET_2015_12_16/benchmarks/rexi_tests/2015_12_15_error_plot_fd/run_rexi_fd_par_m0008_t024_n0032_r0001_a1.txt
#mpiexec.hydra -genv OMP_NUM_THREADS 24 -envall  -np 1 ./build/rexi_fd_par_m_tno_a1  -s 5 -f 1  -g 1 -H 1 -X 1 -Y 1 --compute-error 1 -t 50 -R 4 -C 0.3 -N 32 -U 0 -S 0 --use-fd-for-complex-array 1 --rexi-h 0.8 --timestepping-mode 1 --staggering 0 --rexi-m=8 -C -5.0 #> /gpfs/stfc/local/HCP011/smp03/mxs64-smp03/workspace/SWEET_2015_12_16/benchmarks/rexi_tests/2015_12_15_error_plot_fd/run_rexi_fd_par_m0008_t024_n0032_r0001_a1.txt

