#! /bin/bash


scons --compile-program=swe --gui=enable
scons --compile-program=swe_staggered --gui=enable
scons --compile-program=advection --gui=enable

scons --compile-program=swe --gui=disable
scons --compile-program=swe_staggered --gui=disable
scons --compile-program=advection --gui=disable
