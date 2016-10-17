#! /bin/bash


cd ../

scons --compile-program=swe_staggered --gui=enable --plane-spectral-space=disable

H0=100
CFL=0.1
SIM_TIME=5
SIM_TIMESTEPS=-1

for i in `seq 1 16`; do
	N=$((2**i))
	echo -en "$N:\t"
	./build/example_swe_staggered_spectraldisable -C $CFL -n $N -s 2 -f 0.1 -t $SIM_TIME -T $SIM_TIMESTEPS -l 100 -v 1 -H $H0 | tail -n 1
done
