#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="libpfasst"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libpfasst.a"
PKG_URL_SRC="libpfasst_sweet_2021_11_17.tar.bz2"

config_setup

config_package $@

if [ -z "$MULE_MPIF90" ]; then
	echo_error_exit "MULE_MPIF90 environment variable must be set for libpfasst compilation, see platform configuration"
fi

# Add an empty line to avoid other problems
echo "" >> Makefile.local

# Use SWEET's default Fortran compiler
echo "FC = ${MULE_MPIF90}" >> Makefile.local


#
# Check if gfortran supports -fallow-argument-mismatch and enable it per default
#

TMPDIR="$(mktemp -d)"
echo "" > "$TMPDIR/dummy.f90"
$MULE_MPIF90 -c -fallow-argument-mismatch "$TMPDIR/dummy.f90" -o "$TMPDIR/dummy.o" 2> /dev/null
if [[ $? -eq 0 ]]; then
	echo "$MULE_MPIF90 seems to support -fallow-argument-mismatch, using this per default"
	echo "FFLAGS += -fallow-argument-mismatch" >> Makefile.local
fi
#Comment the next line when installing on linux cluster
rm -rf "${TMPDIR}"


# Disable LTO since this doesn't work on all platforms
echo "LDFLAGS += -fno-lto" >> Makefile.local

# Activate FFTW
echo "USE_FFTW = TRUE" >> Makefile.local

# Activate Verbose make
echo "MKVERBOSE = TRUE" >> Makefile.local

INCDIR="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
mkdir -p "${INCDIR}"

INCDIR_LIBPFASST="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/libpfasst/"
mkdir -p "${INCDIR_LIBPFASST}"

# Add SWEET's include directory
echo "FFLAGS += -I$INCDIR -I$INCDIR_LIBPFASST" >> Makefile.local

# Add LDFLAGS for FFTW
echo "LDFLAGS += -I$INCDIR -I$INCDIR_LIBPFASST" >> Makefile.local

# Set to true to get Debug version
if true; then
	TMPFILE=$(mktemp)
	cp Makefile.local "${TMPFILE}"
	echo "DEBUG = TRUE" > Makefile.local
	echo "" >> Makefile.local
	cat "${TMPFILE}" >> Makefile.local
fi


echo_info "Executing 'make clean'..."
config_exec make clean

config_exec make

echo_info "Installing..."

# Copy modules and change this path if necessary
echo_info cp -v -f ./include/*mod "$INCDIR_LIBPFASST"
cp -v -f ./include/*mod "$INCDIR_LIBPFASST" || echo_error_exit "Failed to install .mod files"

# Copy static library
mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
echo_info cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/" || echo_error_exit "Failed to install libpfasst.a files"


config_success