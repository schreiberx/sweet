#! /bin/bash

BASEDIR=`pwd`


# simtime: 7
# output: 14 times
# ts: about 0.01

RES=16
TS_SIZE=0.0001

STDEXEC="../../../build/swe_rexi_planespectral_omp_libfft_gnu_release -N $RES -M $RES --nonlinear=0 -X 1 -Y 1 -S 1 --staggering=0"
STDEXEC2="../../../build/swe_rexi_planespectral_planedealiasing_omp_libfft_gnu_release -N $RES --nonlinear=0 -X 1 -Y 1 -S 1 --staggering=0"

# normal mode generation method
#for NG in 3 2 1; do
for NG in 1 2 3; do

	# time stepping method
	#for TSM in 2 3; do
	for TSM in 1; do

		#for f in 0.00014584; do
		for f in 1; do #0.1 1 10; do
		for g in 1; do #0.1 1 10; do
		for h in 1; do #0.1 1 10; do
		for r in 1; do #0.1 1 10; do

		# dealiasing
		for d in 0 1; do

				F=$(printf "%02.2f" $f)
				G=$(printf "%02.2f" $g)
				H=$(printf "%02.2f" $h)
				R=$(printf "%02.2f" $r)

				OUTDIR="$BASEDIR/output_ng$NG""_tsm$TSM""_f$F""_g$G""_h$H""_r$R""_dealiasing$d"
				echo "OUTDIR: $OUTDIR"
				mkdir -p "$OUTDIR"

				cd "$OUTDIR"

				if [ "x$d" = "x0" ]; then
					EXEC="$STDEXEC"
				else
					EXEC="$STDEXEC2"
				fi

				EXEC="$EXEC --normal-mode-analysis-generation=$NG"

				EXEC="$EXEC -f $f"
				EXEC="$EXEC -g $g"
				EXEC="$EXEC -H $h"
				EXEC="$EXEC -a $r"

				#EXEC="$EXEC -s 11"
				EXEC="$EXEC -s 0"

				EXEC="$EXEC --timestepping-method=$TSM"
				EXEC="$EXEC --timestepping-order=1"

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


