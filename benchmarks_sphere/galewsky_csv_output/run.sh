#! /bin/bash

BASEDIR="`pwd`"
rm -f ./prog_*

SPHROOT="../../../"
cd "$SPHROOT"

#make clean
scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release || exit 1

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
EXEC="$(ls -1 $SPHROOT/build/swe_sphere_*_release)"
echo "******************************************************************"
echo "$EXEC $PARAMS"
echo "******************************************************************"

$EXEC $PARAMS || exit 1

