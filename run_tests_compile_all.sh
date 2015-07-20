#! /bin/bash

scons --compile-program=advection --gui=disable || exit
scons --compile-program=polvani --gui=disable --spectral-space=enable || exit

scons --compile-program=spectral_visualization --spectral-space=enable --gui=enable || exit

scons --compile-program=swe_nonstaggered_advective --gui=enable || exit
scons --compile-program=swe_nonstaggered_advective --spectral-space=enable || exit
scons --compile-program=swe_nonstaggered_vector_invariant --gui=enable || exit
scons --compile-program=swe_staggered_vector_invariant --gui=enable || exit
#scons --compile-program=swe --gui=enable || exit

scons --compile-program=advection --gui=enable || exit

#scons --compile-program=test_fft --gui=disable --spectral-space=enable --spectral-dealiasing=enable || exit
scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable || exit
