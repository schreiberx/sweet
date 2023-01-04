
# Install binutils since they are outdated on Travis

# prereq: texinfo package
PKGS+=("install_make.sh")
# prereq: make
PKGS+=("install_binutils.sh")

PKGS+=("install_openssl.sh")
PKGS+=("install_cacerts.sh")
# prereq: openssl, cacerts

PKGS+=("install_miniconda.sh")

PKGS+=("install_fftw3.sh")

# prereq: binutils, fftw, python
PKGS+=("install_shtns.sh")
PKGS+=("install_shtns_python.sh")


PKGS+=("install_lapack.sh")

PKGS+=("install_eigen3.sh")

PKGS+=("install_likwid.sh")
PKGS+=("install_numactl.sh")

PKGS+=("install_scons.sh")


PKGS+=("install_mpich.sh")
PKGS+=("install_libpfasst.sh")
