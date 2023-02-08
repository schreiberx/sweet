#! /bin/bash

source ./install_helpers.sh ""

VERSION=4.0.3

PKG_NAME="MPICH"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/mpicc"
#PKG_URL_SRC="mpich-3.3.2.tar.gz"
PKG_URL_SRC="mpich-${VERSION}.tar.gz"


export FC=$F90
export FCFLAGS=$F90FLAGS

unset F90
unset F90FLAGS


#
# Check if gfortran supports -fallow-argument-mismatch and enable it per default
#

if [[ -z "$FC" ]]; then
	# Use gfortran as default compiler
	FC=gfortran
fi

TMPDIR="$(mktemp -d)"
echo "" > "$TMPDIR/dummy.f90"

$FC -c -fallow-argument-mismatch "$TMPDIR/dummy.f90" -o "$TMPDIR/dummy.o" 2> /dev/null
if [[ $? -eq 0 ]]; then
	echo "$FC seems to support -fallow-argument-mismatch, using this per default"
	export FFLAGS="-fallow-argument-mismatch $FFLAGS"
	export FCFLAGS="-fallow-argument-mismatch $FCLAGS"
fi
rm -rf "${TMPDIR}"

config_setup

config_package $@

config_configure	\
	--includedir="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/mpich-${VERSION}"	\
	--enable-shared \
	--enable-cxx \
	--enable-fast=O2 \

config_make_default_install

config_success
