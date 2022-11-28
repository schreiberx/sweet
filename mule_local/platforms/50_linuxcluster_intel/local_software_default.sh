
#MULE_SYSTEM_PACKAGES="libxft-dev texinfo"
MULE_SYSTEM_PACKAGES=""
if [[ "x$DISPLAY" = "x:0" ]]; then
	export MULE_SYSTEM_PACKAGES="$MULE_SYSTEM_PACKAGES libgl1-mesa-dev"
fi

for i in $MULE_SYSTEM_PACKAGES; do
	dpkg -s "$i" >/dev/null 2>&1
	if [[ "$?" != "0" ]]; then
		echo_error "Debian-based system detected and packages missing, please use"
		echo_error "    sudo apt-get install $MULE_SYSTEM_PACKAGES"
		return
	fi
done

#
# Martin Schreiber's laptop (martinium)
#

#PKGS+=("install_gcc7.1.sh")
PKGS+=("install_fftw3.sh")
#PKGS+=("install_eigen3.sh")
PKGS+=("install_lapack.sh")
#PKGS+=("install_openssl.sh")

PKGS+=("install_miniconda.sh")

PKGS+=("install_scons3.sh")
PKGS+=("install_shtns.sh")
PKGS+=("install_shtns_python.sh")

#PKGS+=("install_sdl2.sh")

