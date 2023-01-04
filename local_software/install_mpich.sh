#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="MPICH"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/mpicc"
PKG_URL_SRC="mpich-3.3.2.tar.gz"


export FC=$F90
export FCFLAGS=$F90FLAGS

unset F90
unset F90FLAGS

config_setup

config_package $@


#
# Check if gfortran supports -fallow-argument-mismatch and enable it per default
#
TMPDIR="$(mktemp -d)"
echo "" > "$TMPDIR/dummy.f90"
$FC -c -fallow-argument-mismatch "$TMPDIR/dummy.f90" -o "$TMPDIR/dummy.o" 2> /dev/null
if [[ $? -eq 0 ]]; then
	echo "$FC seems to support -fallow-argument-mismatch, using this per default"
	export FFLAGS="-fallow-argument-mismatch $FFLAGS"
fi
rm -rf "${TMPDIR}"

config_configure

config_make_default_install

config_success
