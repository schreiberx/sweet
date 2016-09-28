#! /bin/bash

for i in 0 1 2 3; do
	scons --program=advection --numa-block-allocator=$i --threading=omp
done

scons --program=swe_rexi --gui=disable --spectral-space=enable || exit

scons --program=advection --gui=disable || exit
#scons --program=polvani --gui=disable --spectral-space=enable || exit

scons --program=spectral_visualization --spectral-space=enable --gui=enable || exit

scons --program=swe_nonstaggered_advective --gui=enable || exit
scons --program=swe_nonstaggered_advective --spectral-space=enable || exit
scons --program=swe_nonstaggered_vector_invariant --gui=enable || exit
scons --program=swe_staggered_vector_invariant --gui=enable || exit

