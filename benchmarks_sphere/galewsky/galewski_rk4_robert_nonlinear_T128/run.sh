#! /bin/bash

BASEDIR="`pwd`"
rm -f ./prog_*

SPHROOT="../../../"
cd "$SPHROOT"

make clean
scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release || exit 1

cd "$BASEDIR"

# h0=g=f=1


OTS=$((120*20))

RES=128

BENCH="galewsky"

#VISCOSITY=100000
VISCOSITY=0

SIMTIME=720000


EXEC="$SPHROOT/build/swe_sphere_*_release -M $RES -C -$TS -o $OTS -u $VISCOSITY -t $SIMTIME --benchmark $BENCH --timestepping-method=ln_erk --timestepping-order=4"
echo "******************************************************************"
echo "$EXEC"
echo "******************************************************************"

$EXEC || exit 1

