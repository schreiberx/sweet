#! /bin/bash

cd ../../
scons --program=swe_sph_and_rexi --mode=release --plane-spectral-space=disable --sphere-spectral-space=enable
#scons --program=swe_sph_and_rexi --mode=debug --plane-spectral-space=disable --sphere-spectral-space=enable
