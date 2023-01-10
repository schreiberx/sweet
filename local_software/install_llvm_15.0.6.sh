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

ln -s clang-${VERSION} "$MULE_SOFTWARE_ROOT/local_software/local/bin/clang++-${VERSION}" || exit 1

# test suits
#config_exec make check-clang check-flang

# Skip tests for debugging
#config_exec make check-clang


if true; then
	#
	# Special hack for quadmath.h
	# We'll simply copy it over from GCC
	#

	# First, determine the GCC compiler (assuming we compile it with gcc)
	LLVM_CXX_COMPILER=$CXX
	if [[ -z "$LLVM_CXX_COMPILER" ]]; then
		LLVM_CXX_COMPILER=g++
	fi
	LLVM_CXX_COMPILER_VERSION=$($LLVM_CXX_COMPILER -dumpfullversion | sed "s/\..*//")

	if [[ -z "$LLVM_CXX_COMPILER_VERSION" ]]; then
		echo_error_exit "Failed to determine GCC version"
	fi

	echo_info "Assuming compiler '${LLVM_CXX_COMPILER}' with version '${LLVM_CXX_COMPILER_VERSION}'"

	if [[ $LLVM_CXX_COMPILER == g++* ]]; then
		echo_info "quadmath.h workaround which is not shipped with LLVM."
		echo_info "We simply copy it to our default include path."

		if [[ ! -z "$LLVM_CXX_COMPILER_VERSION" ]]; then
			cp "/usr/lib/gcc/x86_64-linux-gnu/12/include/quadmath.h" "local/include/" || exit 1
		fi
	fi
fi


config_success

