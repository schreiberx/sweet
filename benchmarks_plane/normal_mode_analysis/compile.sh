#! /bin/bash

cd ../../
scons --program=swe_rexi --mode=release --plane-spectral-space=enable --sphere-spectral-space=disable --plane-spectral-dealiasing=disable
scons --program=swe_rexi --mode=release --plane-spectral-space=enable --sphere-spectral-space=disable --plane-spectral-dealiasing=enable
