#! /bin/bash


export SWEET_ADD_COMPILE_OPTIONS="--mkl=enable --compiler=intel"

for i in linear_gaussian_*; do

	cd "$i"
	./run.sh
	cd ".."

done
