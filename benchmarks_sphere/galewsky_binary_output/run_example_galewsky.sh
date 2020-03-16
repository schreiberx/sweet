#! /bin/bash


#Reproduce Results from Fig. 4 in
#Galewsky, J., Scott, R. K., & Polvani, L. M. (2004). An initial-value problem for testing numerical models of the global shallow-water equations. Tellus, Series A: Dynamic Meteorology and Oceanography, 56(5), 429–440. https://doi.org/10.1111/j.1600-0870.2004.00071.x


#
# 1) Compile with
#
#scons --program=swe_sphere --gui=enable --sphere-spectral-space=enable --compiler=gnu --threading=omp --mode=release
P=`pwd`
cd ../../
scons --program=swe_sphere --sphere-spectral-space=enable --compiler=gnu --threading=omp --mode=release


#
# 2) Run with
#
cd "$P"
../../build/swe_sphere_COMP_plspec_pldeal_spspec_spdeal_numa2_fft_gnu_thomp_release -M 342 --benchmark-name=galewsky --dt=30 --timestepping-method=l_cn_n_erk --timestepping-order=2 --semi-lagrangian-approximate-sphere-geometry=1 -t $((8*24*60*60))  -o $((60*60*6)) --output-file-mode=bin


