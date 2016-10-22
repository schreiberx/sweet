#! /bin/bash

BASEDIR="`pwd`"
rm -f ./prog_*

SPHROOT="../../../"
cd "$SPHROOT"

make clean
scons --program=swe_sph_and_rexi --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release || exit 1

cd "$BASEDIR"

# h0=g=f=1

# Stable time step size for T64
TS=120

# unstable for explicit RK4 time stepping
#TS=0.01

OTS=$((120*20))

RES=64
#RES=16

USE_REXI=0
REXI_M=0	# Deactivate REXI
REXI_H=0.2
REXI_HALF_POLES=1
REXI_EXTENDED_MODES=4

BENCH_ID=1
USE_ROBERT_FUNS=1

NONLINEAR=1
VISCOSITY=100000

SIMTIME=720000


$SPHROOT/build/swe_sph_and_rexi_spherespectral_spheredealiasing_*_release -M $RES --nonlinear $NONLINEAR -C -$TS -o $OTS -u $VISCOSITY -t $SIMTIME --rexi=$USE_REXI --rexi-m=$REXI_M --rexi-h=$REXI_H --rexi-half=$REXI_HALF_POLES -s $BENCH_ID --use-robert-functions $USE_ROBERT_FUNS --rexi-ext-modes $REXI_EXTENDED_MODES || exit 1


#mv prog_* "$BASEDIR"

$SPHROOT/plot_csv.py prog_h_*.csv
$SPHROOT/create_mp4.sh prog_h out_prog_h.mp4

$SPHROOT/plot_csv.py prog_eta_*.csv
$SPHROOT/create_mp4.sh prog_eta out_prog_eta.mp4
