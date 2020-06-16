#! /bin/bash

source ./install_helpers.sh ""

# Setup environment in case it has been missed
#source "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate"

PKG_NAME="SHTNS_python"

PYTHONVERSION=$(python3 -c "import sys;print(str(sys.version_info.major)+\".\"+str(sys.version_info.minor),end='')")
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/python$PYTHONVERSION/site-packages/shtns.py"

PKG_URL_SRC="shtns-3.4.tar.gz"

config_setup

config_package $@

if [ "#$TRAVIS" != "#" ]; then
	echo_info "Detected Travis"
	echo_info "Disabling SIMD because of Travis"
	CONFIGURE_EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-simd"
fi

CONFIGURE_EXTRA_FLAGS+=" --enable-ishioka"


echo_info_hline
echo_info "SHTNS Python noOpenMP:"
# Python, no OpenMP

# Special flag for sk2 (@ CAPS hardware)
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
       export CFLAGS="$CFLAGS -march=skylake"
fi

config_configure --enable-python --disable-openmp $CONFIGURE_EXTRA_FLAGS

# Special flag for sk2 (@ CAPS hardware)
pwd
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
	sed -i "s/-march=native/-march=skylake/" "Makefile"
fi

config_make_clean
config_make_default
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

echo_info_hline
echo_info "SHTNS Python OpenMP:"
# Python, OpenMP
config_configure --enable-python --enable-openmp $CONFIGURE_EXTRA_FLAGS

# Special flag for sk2 (@ CAPS hardware)
pwd
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
	sed -i "s/-march=native/-march=skylake/" "Makefile"
fi

config_make_clean
config_make_default
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

config_success
