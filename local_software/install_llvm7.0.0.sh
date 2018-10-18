#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="llvm_7.0.0"
PKG_URL_SRC="http://releases.llvm.org/7.0.0/llvm-7.0.0.src.tar.xz"

config_setup

config_package $@

# Clang
cd llvm/tools
config_download_and_extract http://releases.llvm.org/7.0.0/cfe-7.0.0.src.tar.xz
cd ../../

# Clang tools
cd llvm/tools/clang/tools
config_download_and_extract http://releases.llvm.org/7.0.0/clang-tools-extra-7.0.0.src.tar.xz
cd ../../../..

# Compiler-RT
cd llvm/projects
config_download_and_extract http://releases.llvm.org/7.0.0/compiler-rt-7.0.0.src.tar.xz
cd ../..

cd llvm/tools
config_download_and_extract https://github.com/flang-compiler/flang/archive/flang_20180921.tar.gz
cd ../../

# Compile
mkdir -p build
cd build
config_exec cmake -G "Unix Makefiles" ../llvm
config_exec make

# test suits
config_exec make check-clang

#for i in g++ gcc gcc-ar gcc-nm gcc-ranlib gfortran gcov gcov-tool gfortran; do
#	ln -sf "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/$i-8.2" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/$i" || echo_error_exit "FAILED to create sym links"
#done

config_success

