#! /bin/bash

BASEDIR="`pwd`"
rm -f ./prog_*

SWEETROOT="../../../"
cd "$SWEETROOT"

make clean
scons --program=swe_sph_and_rexi --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release $SWEET_ADD_COMPILE_OPTIONS

cd "$BASEDIR"

# h0=g=f=1

# Stable time step size for T64
TS=0.001

# Instable with explicit methods
TS=0.01

# unstable for explicit RK4 time stepping
#TS=0.01


RES=64
#RES=16

OTS=0.01

USE_REXI=1
REXI_M=$((256))
REXI_H=0.2
REXI_HALF_POLES=1
REXI_EXTENDED_MODES=4

G=1	# gravity
H=1	# avg height
F=0	# omega
R=1	# radius

BENCH_ID=3
USE_ROBERT_FUNS=0

NONLINEAR=0
VISCOSITY=0

SIMTIME=4

EXEC="$SWEETROOT/build/swe_sph_and_rexi_*_release -g $G -H $H -f $F -a $R -M $RES $NONLINEAR -C -$TS -o $OTS -u $VISCOSITY -t $SIMTIME --nonlinear=$NONLINEAR --rexi=$USE_REXI --rexi-m=$REXI_M --rexi-h=$REXI_H --rexi-half=$REXI_HALF_POLES -s $BENCH_ID --use-robert-functions $USE_ROBERT_FUNS --rexi-ext-modes $REXI_EXTENDED_MODES"
$EXEC


$SWEETROOT/plot_csv.py prog_h_*.csv
$SWEETROOT/create_mp4.sh prog_h out_prog_h.mp4
