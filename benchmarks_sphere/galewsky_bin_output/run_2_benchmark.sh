#! /bin/bash

BASEDIR="`pwd`"
SWEETROOT="../../"

cd "$BASEDIR"

# h0=g=f=1


TS=$((120))
OTS=$((TS*20))

RES=128

BENCH="galewsky"

#VISCOSITY=100000
VISCOSITY=0

SIMTIME=720000


PARAMS=""
#PARAMS+=" --output-file-mode bin"
PARAMS+=" -M $RES"
PARAMS+=" --dt $TS"
PARAMS+=" -o $OTS -u $VISCOSITY -t $SIMTIME --benchmark-name $BENCH --timestepping-method=ln_erk --timestepping-order=4"
PARAMS+=" --output-file-mode=bin"
EXEC="$(ls -1 $SWEETROOT/build/swe_sphere_*_release)"
echo "******************************************************************"
echo "$EXEC $PARAMS"
echo "******************************************************************"

$EXEC $PARAMS || exit 1

