#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo_info_hline
echo_info "PFASST"
echo_info_hline


if [ -e "$MULE_LOCAL_ROOT/../local_software/local/lib/libpfasst.a" ]; then 

	SCONS="scons --program=libpfasst_swe_sphere_expl_sdc --quadmath=disable --libpfasst=enable --sweet-mpi=enable --libsph=enable --numa-block-allocator=0 --plane-spectral-space=disable --sphere-spectral-space=enable --threading=off --libfft=enable --sphere-spectral-dealiasing=enable"
	echo "$SCONS"
	SCONS="scons --program=libpfasst_swe_sphere --quadmath=disable --libpfasst=enable --sweet-mpi=enable --libsph=enable --numa-block-allocator=0 --plane-spectral-space=disable --sphere-spectral-space=enable --threading=off --libfft=enable --sphere-spectral-dealiasing=enable"
	echo "$SCONS"
	$SCONS || exit

fi
