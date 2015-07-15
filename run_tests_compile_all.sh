#! /bin/bash

scons --compile-program=advection --gui=disable || exit
scons --compile-program=polvani --gui=disable --spectral-space=enable || exit

scons --compile-program=spectral_visualization --spectral-space=enable --gui=enable || exit

scons --compile-program=swe_nonstaggered_covariant --gui=enable || exit
scons --compile-program=swe_staggered_covariant --gui=enable || exit
scons --compile-program=swe --gui=enable || exit

scons --compile-program=advection --gui=enable || exit

scons --compile-program=swe --gui=enable || exit
scons --compile-program=swe --gui=disable || exit

scons --compile-program=test_fft --gui=disable --spectral-space=enable --spectral-dealiasing=enable || exit
scons --compile-program=test_spectral_ops --gui=disable --spectral-space=enable || exit
