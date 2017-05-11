#! /bin/bash

BASEDIR=`pwd`


# simtime: 7
# output: 14 times
# ts: about 0.01

# simtime: 1.5 days (129600 seconds)
# output: 12 times (10800 seconds)
# ts: about 200

RESO=32
RES=$RESO
TS_SIZE=0.000000001
NUM_TS=10

TS_SIZE=1
NUM_TS=1

FSPHERE=1

STDEXEC="../../../build/swe_sph_and_rexi_spherespectral_spheredealiasing_omp_libfft_gnu_release --nonlinear=0 --rexi-m=16 -F $FSPHERE"
#STDEXEC="../../../build/swe_sph_and_rexi_spherespectral_spheredealiasing_omp_libfft_gnu_debug --nonlinear=0 --rexi-m=16"

# normal mode generation method
#for NG in 1 3; do
for NG in 3; do

	# reduce resolution for physical space
	if [ "x$NG" == "x1" ]; then
		RES=$(($RESO/2))
	fi

	#
	# Use solver which would include Coriolis effect or not
	#
#	for COREFF in 0 1; do
	for COREFF in 1; do

#		for EXT_MODES in 0 2 4; do
		for EXT_MODES in 0; do

			# time stepping method
			#for TSM in 2 3; do
			# 1: Explicit RK
			# 2: Leapfrog
			# 3: Implicit Euler (order 1) and Crank-Nicolson (order 2)
			# 100: REXI
			#for TSM in 1 2 3; do
			#for TSM in 1 100; do
			for TSM in 1; do # 100; do

				# time stepping order
				#for TSMO in 1 2; do
				for TSMO in 1; do

					# pde id to analyse
					#for PDEID in 0 1 2; do
					#for PDEID in 0 1 2; do
					#for PDEID in 0 1 2; do
					for PDEID in 1; do

						# robert formulation?
						#for ROB in 0 1; do
						for ROB in 1; do


							if [ $TSM -ge 100 ]; then
								REXIM_SET="128 256 512"
							else
								REXIM_SET="1"
							fi

							for REXI_M in $REXIM_SET; do

								#for f in 0 1; do #0.1 1 10; do
								for f in 0.00014584 ; do #0.1 1 10; do
								for g in 1; do #0.1 1 10; do
								for h in 100000; do #0.1 1 10; do
								for r in 6371220; do #0.1 1 10; do

									F=$(printf "%05.5f" $f)
									G=$(printf "%05.5f" $g)
									H=$(printf "%05.5f" $h)
									R=$(printf "%05.5f" $r)

									OUTDIR="$BASEDIR/output_pde$PDEID""_ts$NUM_TS""_rob$ROB""_numModeGen$NG""_timeStMethod$TSM""_tsOrder$TSMO""_rexim$REXI_M""_f$F""_grav$G""_height$H""_radius$R""_extmodes$EXT_MODES""_coreff$COREFF"

									echo "OUTDIR: $OUTDIR"
									mkdir -p "$OUTDIR"

									cd "$OUTDIR"

									EXEC="$STDEXEC --normal-mode-analysis-generation=$NG --rexi-ext-modes=$EXT_MODES --use-coriolis-formulation=$COREFF --rexi-m=$REXI_M"

									EXEC="$EXEC -f $f"
									EXEC="$EXEC -g $g"
									EXEC="$EXEC -H $h"
									EXEC="$EXEC -a $r"

									#EXEC="$EXEC -s 11"
									EXEC="$EXEC -s 0"

									EXEC="$EXEC --timestepping-method=$TSM"
									EXEC="$EXEC --timestepping-order=$TSMO"

									EXEC="$EXEC --pde-id=$PDEID"
									EXEC="$EXEC --use-robert-functions=$ROB"

									EXEC="$EXEC -M $RES"
									EXEC="$EXEC -C -$TS_SIZE"

									echo "$EXEC"

									$EXEC

									cd "$BASEDIR"
								done
								done
								done
								done
							done
						done
					done
				done
			done
		done
	done
done


echo "*********************************************"
echo "* TS_SIZE: $TS_SIZE"
echo "* NUM_TS: $NUM_TS"
echo "*********************************************"
