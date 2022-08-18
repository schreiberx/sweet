#!/bin/bash

## How to use
## ./postprocessing_consolidate.py iteration/timestep physical/spectral var1 name_var1 vals_var1 var2 name_var2 vals_var2 ...
## iteration/timestep: variable for x axis
## physical/spectral: error computation
## varX name_varX vals_varX: only plot simulations with name_varX = vals_varX; no filter is applied to non listed parameters

## plot errors along iterations for given numbers of coarse slices
./postprocessing_consolidate.py iteration physical var1 runtime.parareal_coarse_slices 5 var2 runtime.parareal_coarse_timestep_size 15 60 120
./postprocessing_consolidate.py iteration physical var1 runtime.parareal_coarse_slices 10 var2 runtime.parareal_coarse_timestep_size 15 60 120

## plot errors at fixed iteration in function of the coarse timestep size for given numbers of coarse slices
./postprocessing_consolidate.py timestep physical var1 runtime.parareal_coarse_slices 5 var2 runtime.iteration 1
./postprocessing_consolidate.py timestep physical var1 runtime.parareal_coarse_slices 10 var2 runtime.iteration 1
