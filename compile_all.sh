#! /bin/bash


scons --compile-program=swe --gui=enable
scons --compile-program=swe_staggered --gui=enable
scons --compile-program=advection --gui=enable
