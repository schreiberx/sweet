#! /bin/bash

dir=`pwd`
echo $dir
cd $dir

source install_autoconf.sh
cd $dir
#source install_automake.sh # zip corrupted
cd $dir
#source install_eigen3.sh # zip corrupted
cd $dir
source install_fftw3.sh
cd $dir
#source install_gcc5.3.sh
cd $dir
#source install_gcc7.1.sh
cd $dir
source install_lapack.sh
cd $dir
source install_libfreetype.sh
cd $dir
#source install_libpfasst.sh #needs password
cd $dir
#source install_libpng.sh
cd $dir
source install_likwid.sh
cd $dir
source install_numa.sh
cd $dir
source install_openmpi.sh
cd $dir
source install_pfasst++.sh
cd $dir
source install_python2.sh
cd $dir
source install_python_matplotlib.sh
cd $dir
source install_python_mpmath.sh
cd $dir
source install_python_pip.sh
cd $dir
source install_python_scipy.sh
cd $dir
source install_rdic.sh
cd $dir
source install_scons.sh
cd $dir
source install_sdl2.sh
cd $dir
source install_shtns.sh
cd $dir
