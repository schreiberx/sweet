#! /bin/bash

PARAMS=""
PARAMS+=" --tool=memcheck"
PARAMS+=" --leak-check=full --show-reachable=yes --track-origins=yes --trace-children=yes"
#PARAMS+=" --tool=memcheck"
#PARAMS+=" --suppressions=mpich-valgrind.supp"
#PARAMS+=" --suppressions=/home/schreibm/local/openmpi-1.6.3/share/openmpi/openmpi-valgrind.supp"
echo valgrind $PARAMS $@
valgrind $PARAMS $@

