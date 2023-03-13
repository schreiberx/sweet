#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"


for PLANE_SPECTRAL_SPACE in enable disable; do
	for REXI_THREAD in enable disable; do
		for SWEET_MPI in enable disable; do
			for THREADING in off omp; do
				echo $REXI_THREAD

				echo
				SCONS="scons --program=swe_sphere --gui=disable --sphere-spectral-space=enable --mode=debug --gui=enable"
			      	SCONS+=" --plane-spectral-space=$PLANE_SPECTRAL_SPACE"
			       	SCONS+=" --rexi-thread-parallel-sum=$REXI_THREAD"
			       	SCONS+=" --sweet-mpi=$SWEET_MPI"
			       	SCONS+=" --threading=$THREADING"
				echo "$SCONS"
				$SCONS  || exit

				mule.benchmark.cleanup_all || exit 1
			done
		done
	done
done
