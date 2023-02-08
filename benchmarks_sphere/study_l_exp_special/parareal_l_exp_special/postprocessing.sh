#!/bin/bash

## How to use
## ./postprocessing_consolidate.py iteration/timestep physical/spectral var1 name_var1 vals_var1 var2 name_var2 vals_var2 ...
## iteration/timestep: variable for x axis
## physical/spectral: error computation
## varX name_varX vals_varX: only plot simulations with name_varX = vals_varX; no filter is applied to non listed parameters

plot_timings=0

for timeslice_size in {10,5}; do
    for space in {physical,spectral}; do

        ## plot errors along iterations for given numbers of coarse slices and given coarse timestep sizes
        ./postprocessing_consolidate.py iteration $space $plot_timings var1 runtime.parareal_coarse_slices $timeslice_size var2 runtime.parareal_coarse_timestep_size 30 60 120

        ## plot errors at fixed iteration in function of the coarse timestep size for given numbers of coarse slices
        for niter in {1,2,5,7,10}; do
            if [ $niter -gt $timeslice_size ]; then
                continue;
            fi;
            ./postprocessing_consolidate.py timestep $space $plot_timings var1 runtime.parareal_coarse_slices $timeslice_size var2 runtime.iteration $niter
        done
    done
done

