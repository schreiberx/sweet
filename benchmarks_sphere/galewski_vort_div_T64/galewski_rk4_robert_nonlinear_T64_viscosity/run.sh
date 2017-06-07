#! /bin/bash

BASEDIR="`pwd`"
rm -f ./prog_*

SPHROOT="../../../"
cd "$SPHROOT"

make clean
scons --program=swe_sphere_rexi --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release || exit 1

cd "$BASEDIR"

# h0=g=f=1

# Stable time step size for T64
TS=120

# unstable for explicit RK4 time stepping
#TS=0.01

OTS=$((120*20))

RES=64
#RES=16

REXI_M=128
REXI_H=0.2
REXI_HALF_POLES=1
REXI_EXTENDED_MODES=4

BENCH_ID=1

NONLINEAR=1
VISCOSITY=100000

SIMTIME=720000

EXEC="$SPHROOT/build/swe_sphere_rexi_spherespectral_spheredealiasing_*_release -M $RES -C -$TS -o $OTS -u $VISCOSITY -t $SIMTIME --rexi-m=$REXI_M --rexi-h=$REXI_H --rexi-half=$REXI_HALF_POLES -s $BENCH_ID --rexi-ext-modes $REXI_EXTENDED_MODES --timestepping-method=ln_erk --timestepping-order=4"
echo "******************************************************************"
echo "$EXEC"
echo "******************************************************************"

$EXEC || exit 1

