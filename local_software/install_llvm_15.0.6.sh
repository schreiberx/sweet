#! /bin/bash

source ./install_helpers.sh ""

VERSION=15
VERSION_FULL=15.0.6

# Setup environment in case it has been missed (e.g. Conda was installed right before)
if [[ -e "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate" ]]; then
	source "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate" || exit 1
fi

PKG_NAME="llvm-${VERSION_FULL}"
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/llvmorg-${VERSION_FULL}.tar.gz"
# Package based on GIT release package
#PKG_URL_SRC="https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-15.0.6.tar.gz"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/clang-${VERSION}"


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

echo_info "Creating symlink for 'clang++-${VERSION}'"

ln -s clang "local/bin/clang++-${VERSION}" || exit 1

# test suits
#config_exec make check-clang check-flang

# Skip tests for debugging
#config_exec make check-clang

config_success

