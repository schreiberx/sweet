#! /bin/bash


#export SWEET_ADD_COMPILE_OPTIONS="--mkl=enable --compiler=intel"
export SWEET_ADD_COMPILE_OPTIONS="--threading=off --rexi-parallel-sum=enable"


for i in linear_gaussian_*; do

	cd "$i"
	./run.sh || exit
	cd ".."

done
