#! /bin/bash

DIR=$(dirname $0)

NAME="$1"
$DIR/normal_modes_compute.py  "$NAME"
$DIR/normal_modes_plot.py  "$NAME"_evalues_real.csv  "$NAME"_evalues_imag.csv  && eog output.png
