
# Use this environment variable to add further packages
for i in $SWEET_LOCAL_SOFTWARE_PRE; do
	PKGS+=($i)
done

PKGS+=("install_miniconda.sh")

PKGS+=("install_fftw3.sh")
PKGS+=("install_shtns.sh")
PKGS+=("install_shtns_python.sh")

PKGS+=("install_scons.sh")

PKGS+=("install_mpich.sh")

PKGS+=("install_eigen3.sh")
PKGS+=("install_lapack.sh")

PKGS+=("install_numactl.sh")

#PKGS+=("install_libpfasst.sh")
PKGS+=("install_libpfasst_debug.sh")
PKGS+=("install_xbraid.sh")

PKGS+=("install_sdl2.sh")


# Use this environment variable to add further packages
for i in $SWEET_LOCAL_SOFTWARE_POST; do
	PKGS+=($i)
done
