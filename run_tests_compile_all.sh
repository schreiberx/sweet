#! /bin/bash

for i in 0 1 2 3; do
	scons --program=advection --numa-block-allocator=$i --threading=omp --mode=debug
done

scons --program=burgers --parareal=serial --plane-spectral-space=enable --mode=debug || exit

scons --program=swe_rexi --gui=disable --plane-spectral-space=enable --mode=debug || exit

scons --program=advection --gui=disable --mode=debug || exit
#scons --program=polvani --gui=disable --plane-spectral-space=enable || exit

scons --program=spectral_visualization --plane-spectral-space=enable --gui=enable --mode=debug || exit

scons --program=swe_nonstaggered_advective --gui=enable --mode=debug || exit
scons --program=swe_nonstaggered_advective --plane-spectral-space=enable --mode=debug || exit
scons --program=swe_nonstaggered_vector_invariant --gui=enable --mode=debug || exit
scons --program=swe_staggered_vector_invariant --gui=enable --mode=debug || exit

