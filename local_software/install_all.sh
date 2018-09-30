#! /bin/bash

source ./config.sh ""
source ./env_vars.sh ""

PKGS=()

if [ "${HOSTNAME:0:8}" == "cheyenne" ]; then

	echo "***********************************"
	echo "Cheyenne system detected"
	echo "***********************************"

	#PKGS+=("install_autoconf.sh")
	#PKGS+=("install_automake.sh")
	PKGS+=("install_fftw3.sh")
	PKGS+=("install_eigen3.sh")
	#PKGS+=("install_gcc5.3.sh")
	#PKGS+=("install_gcc7.1.sh")
	PKGS+=("install_lapack.sh")
	#PKGS+=("install_libfreetype.sh")
	#PKGS+=("install_libpfasst.sh")
	#PKGS+=("install_libpng.sh")
	#PKGS+=("install_likwid.sh")
	#PKGS+=("install_numa.sh")
	#PKGS+=("install_openmpi.sh")
	#PKGS+=("install_pfasst++.sh")
	#PKGS+=("install_python3_pip.sh")
	PKGS+=("install_python3.sh")
	#PKGS+=("install_rdic.sh")
	PKGS+=("install_scons3.sh")
	#PKGS+=("install_sdl2.sh")
	PKGS+=("install_shtns.sh")

elif [ "${HOSTNAME:0:5}" == "mpp2-" ]; then
	#
	# LRZ Munich
	# CoolMUC Cluster
	# mpp2 partition
	#
	#PKGS+=("install_gcc7.1.sh")
	PKGS+=("install_fftw3.sh")
	#PKGS+=("install_eigen3.sh")
	PKGS+=("install_lapack.sh")
	PKGS+=("install_openssl.sh")
	PKGS+=("install_python3.sh")
	PKGS+=("install_scons3.sh")
	PKGS+=("install_shtns.sh")
	PKGS+=("install_shtns_python.sh")

elif [ "${HOSTNAME:0:9}" == "martinium" ]; then
	#
	# Martin Schreiber's laptop (martinium)
	#
	#PKGS+=("install_gcc7.1.sh")
	PKGS+=("install_fftw3.sh")
	#PKGS+=("install_eigen3.sh")
	PKGS+=("install_lapack.sh")
	#PKGS+=("install_openssl.sh")
	#PKGS+=("install_python3.sh")
	PKGS+=("install_scons3.sh")
	PKGS+=("install_shtns.sh")
	PKGS+=("install_shtns_python.sh")
#	PKGS+=("install_numa.sh")

elif [ "${HOSTNAME}" == "mac-login-amd" ]; then
	PKGS+=("install_gcc7.1.sh")
	PKGS+=("install_fftw3.sh")
	PKGS+=("install_eigen3.sh")
	PKGS+=("install_lapack.sh")
	PKGS+=("install_python3.sh")
	PKGS+=("install_scons3.sh")
	PKGS+=("install_shtns.sh")

elif [ "${HOSTNAME:0:8}" == "guepardo" ] || [ "${HOSTNAME:0:5}" == "green" ]; then

	echo "***********************************"
	echo "USP system detected at USP-BR"
	echo "***********************************"

	#PKGS+=("install_autoconf.sh")
	#PKGS+=("install_automake.sh")
	PKGS+=("install_fftw3.sh")
	PKGS+=("install_eigen3.sh")
	#PKGS+=("install_gcc5.3.sh")
	PKGS+=("install_gcc7.2.sh")
	PKGS+=("install_lapack.sh")
	#PKGS+=("install_libfreetype.sh")
	#PKGS+=("install_libpfasst.sh")
	PKGS+=("install_libpng.sh")
	PKGS+=("install_likwid.sh")
	PKGS+=("install_numa.sh")
	PKGS+=("install_openmpi.sh")
	#PKGS+=("install_pfasst++.sh")
	PKGS+=("install_python3.sh")
	PKGS+=("install_python3_pip.sh")
	#PKGS+=("install_rdic.sh")
	PKGS+=("install_scons3.sh")
	PKGS+=("install_sdl2.sh")
	PKGS+=("install_shtns.sh")
else
	echo "***********************************"
	echo "System not detected, aborting"
	echo "***********************************"
	exit 1
fi

echo "Installing packages: ${PKGS[@]}"

for D in "${PKGS[@]}"; do
	./$D || exit 1
done
