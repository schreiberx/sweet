#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo_info_hline
echo_info "PFASST-EXPL-SDC"
echo_info_hline


make clean

SCONS="scons --program=libpfasst_swe_sphere_mlsdc --quadmath=enable --libpfasst=enable --sweet-mpi=enable --libsph=enable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=off --libfft=enable --mode=debug"
echo "$SCONS"
$SCONS || exit

NSWEEPS_COARSE_1="./build/libpfasst_swe_sphere_mlsdc_*_debug -M 64 -t 360 --benchmark-name=gaussian_bump --dt=180 --libpfasst-nlevels 2 --libpfasst-niters 1 --libpfasst-nnodes 3 --output-file-mode=bin --libpfasst-nsweeps 1 -o 360"
echo "$NSWEEPS_COARSE_1"
$NSWEEPS_COARSE_1 || exit

NSWEEPS_COARSE_2="./build/libpfasst_swe_sphere_mlsdc_*_debug -M 64 -t 360 --benchmark-name=gaussian_bump --dt=180 --libpfasst-nlevels 2 --libpfasst-niters 1 --libpfasst-nnodes 3 --output-file-mode=bin --libpfasst-nsweeps 2,1 -o 360"
echo "$NSWEEPS_COARSE_2"
$NSWEEPS_COARSE_2 || exit

mule.benchmark.cleanup_all || exit 1

