#! /usr/bin/env bash


./postprocessing_pickle.py || exit 1

./postprocessing_consolidate_prog_phi_pert.py || exit 1
./postprocessing_consolidate_prog_div.py || exit 1
./postprocessing_consolidate_prog_vrt.py || exit 1

