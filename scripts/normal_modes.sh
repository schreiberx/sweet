#! /bin/bash

DIR=$(dirname $0)

NAME="$1"
$DIR/normal_modes_compute.py  "$NAME"
$DIR/normal_modes_plot_and_analyse.py  "$NAME"_evalues_complex.csv
#eog output.png
