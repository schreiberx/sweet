#! /bin/bash

source ./install_helpers.sh "" || exit 1

# Name of package
PKG_NAME="llvm_trunk"

# Path to one file of installed package to test for existing installation
#PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/gcc-8.2"

# URL to source code to fetch it
#PKG_URL_SRC="gcc-8.2.0.tar.xz"

config_package $@ || exit 1

config_exec svn co http://llvm.org/svn/llvm-project/llvm/trunk llvm
cd llvm/tools

#Clang
cd llvm/tools
svn co http://llvm.org/svn/llvm-project/cfe/trunk clang
cd ../..

#Clang tools
cd llvm/tools/clang/tools
svn co http://llvm.org/svn/llvm-project/clang-tools-extra/trunk extra
cd ../../../..

# Compiler-RT
cd llvm/projects
svn co http://llvm.org/svn/llvm-project/compiler-rt/trunk compiler-rt
cd ../..

# Compile
mkdir build (in-tree build is not supported)
cd build
cmake -G "Unix Makefiles" ../llvm
make


#for i in g++ gcc gcc-ar gcc-nm gcc-ranlib gfortran gcov gcov-tool gfortran; do
#	ln -sf "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/$i-8.2" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/$i" || echo_error_exit "FAILED to create sym links"
#done

config_success || exit 1

