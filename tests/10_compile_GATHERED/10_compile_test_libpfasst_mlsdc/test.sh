#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo_info_hline
echo_info "PFASST-MLSDC"
echo_info_hline


if [ -e "$MULE_LOCAL_ROOT/../local_software/local/lib/libpfasst.a" ]; then 

	SCONS="scons --program=libpfasst_swe_sphere_mlsdc --quadmath=disable --libpfasst=enable --sweet-mpi=enable --libsph=enable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=off --libfft=enable --sphere-spectral-dealiasing=enable"
	echo "$SCONS"
	$SCONS || exit

fi
