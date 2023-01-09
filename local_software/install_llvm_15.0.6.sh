#! /bin/bash

source ./install_helpers.sh ""

# Setup environment in case it has been missed (e.g. Conda was installed right before)
if [[ -e "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate" ]]; then
	source "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate" || exit 1
fi

#
# Package based on GIT release package
# https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.6/clang-15.0.6.src.tar.xz
#

PKG_NAME="llvm-15.0.6"
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/llvmorg-15.0.6.tar.gz"
#PKG_URL_SRC="https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-15.0.6.tar.gz"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/clang-15"


config_setup

config_package $@

# Compile
mkdir -p build
cd build

# Project 'flang' is not yet included since 'flang' is not supported by MPICH, yet
config_exec cmake \
	-DLLVM_ENABLE_PROJECTS="clang;openmp"	\
	-DCMAKE_BUILD_TYPE=Release \
	-DLLVM_TARGETS_TO_BUILD=Native	\
	-DLLVM_BUILD_TOOLS=OFF \
	-DLLVM_BUILD_UTILS=OFF  \
	-G "Unix Makefiles" \
	-DCMAKE_INSTALL_PREFIX:PATH=$SWEET_LOCAL_SOFTWARE_DST_DIR \
	../llvm

config_make_default
config_exec make install

# test suits
#config_exec make check-clang check-flang
config_exec make check-clang

config_success

