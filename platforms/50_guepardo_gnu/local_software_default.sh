

SWEET_SYSTEM_PACKAGES="libxft-dev libssl-dev texinfo"
if [[ "$DISPLAY" = ":0" ]]; then
	export SWEET_SYSTEM_PACKAGES="$SWEET_SYSTEM_PACKAGES libgl1-mesa-dev libxext-dev"
fi

for i in $SWEET_SYSTEM_PACKAGES; do
	dpkg -s "$i" >/dev/null 2>&1
	if [[ "$?" != "0" ]]; then
		echo_error "Debian-based system detected and packages missing, please use"
		echo_error "    sudo apt-get install $SWEET_SYSTEM_PACKAGES"
		return 2>/dev/null
		exit 1
	fi
done


#PKGS+=("install_autoconf.sh")
PKGS+=("install_make.sh")
PKGS+=("install_gcc8.2.sh")

PKGS+=("install_cmake.sh")
#PKGS+=("install_automake.sh")

PKGS+=("install_fftw3.sh")
PKGS+=("install_eigen3.sh")
#PKGS+=("install_gcc5.3.sh")
#PKGS+=("install_gcc7.2.sh")
PKGS+=("install_lapack.sh")

PKGS+=("install_likwid.sh")
PKGS+=("install_numactl.sh")
#PKGS+=("install_pfasst++.sh")

PKGS+=("install_python3.sh")
PKGS+=("install_scons3.sh")
PKGS+=("install_shtns.sh")
PKGS+=("install_shtns_python.sh")


PKGS+=("install_openmpi.sh")
PKGS+=("install_libpfasst.sh")

PKGS+=("install_libfreetype.sh")
PKGS+=("install_libpng.sh")
PKGS+=("install_sdl2.sh")
