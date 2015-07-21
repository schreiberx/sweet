#! /bin/bash

scons --program=advection --gui=disable || exit
scons --program=polvani --gui=disable --spectral-space=enable || exit

scons --program=spectral_visualization --spectral-space=enable --gui=enable || exit

scons --program=swe_nonstaggered_advective --gui=enable || exit
scons --program=swe_nonstaggered_advective --spectral-space=enable || exit
scons --program=swe_nonstaggered_vector_invariant --gui=enable || exit
scons --program=swe_staggered_vector_invariant --gui=enable || exit

scons --unit-test=test_advection --gui=enable || exit
scons --unit-test=test_spectral_ops --gui=disable --spectral-space=enable || exit
