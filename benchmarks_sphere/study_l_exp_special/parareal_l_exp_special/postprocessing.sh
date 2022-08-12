#!/bin/bash

## plot errors along iterations for given numbers of coarse slices
./postprocessing_consolidate.py iteration var1 runtime.parareal_coarse_slices 5 var2 runtime.parareal_coarse_timestep_size 15 60 120
./postprocessing_consolidate.py iteration var1 runtime.parareal_coarse_slices 10 var2 runtime.parareal_coarse_timestep_size 15 60 120

## plot errors at fixed iteration in function of the coarse timestep size for given numbers of coarse slices
./postprocessing_consolidate.py timestep var1 runtime.parareal_coarse_slices 5 var2 runtime.iteration 1
./postprocessing_consolidate.py timestep var1 runtime.parareal_coarse_slices 10 var2 runtime.iteration 1
