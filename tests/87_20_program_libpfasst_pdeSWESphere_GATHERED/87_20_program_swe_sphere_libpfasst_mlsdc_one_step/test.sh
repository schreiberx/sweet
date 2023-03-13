#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo_info_hline
echo_info "PFASST-EXPL-SDC"
echo_info_hline


make clean

SCONS="scons --program=programs/libpfasst/pde_sweSphere_mlsdc --libpfasst=enable --sweet-mpi=enable --libsph=enable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=off --libfft=enable --mode=debug"
echo "$SCONS"
$SCONS || exit 1

NSWEEPS_COARSE_1="./build/programs/libpfasst/pde_sweSphere_mlsdc_COMP*_debug -M 64 -t 360 --benchmark-name=gaussian_bump --dt=180 --libpfasst-nlevels 2 --libpfasst-niters 1 --libpfasst-nnodes 3 --output-file-mode=bin --libpfasst-nsweeps 1 -o 360"
echo "$NSWEEPS_COARSE_1"
$NSWEEPS_COARSE_1 || exit 1

NSWEEPS_COARSE_2="./build/programs/libpfasst/pde_sweSphere_mlsdc_COMP*_debug -M 64 -t 360 --benchmark-name=gaussian_bump --dt=180 --libpfasst-nlevels 2 --libpfasst-niters 1 --libpfasst-nnodes 3 --output-file-mode=bin --libpfasst-nsweeps 2,1 -o 360"
echo "$NSWEEPS_COARSE_2"
$NSWEEPS_COARSE_2 || exit 1

mule.benchmark.cleanup_all || exit 1

