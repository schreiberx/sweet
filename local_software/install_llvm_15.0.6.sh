#! /bin/bash

source ./install_helpers.sh ""

#
# Package based on GIT release package
# https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.6/clang-15.0.6.src.tar.xz
#

PKG_NAME="llvm-15.0.6"
#PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/llvmorg-9.0.1.tar.gz"
PKG_URL_SRC="https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-15.0.6.tar.gz"

config_setup

config_package $@

# Compile
mkdir -p build
cd build
config_exec cmake -DLLVM_ENABLE_PROJECTS="clang;flang;polly;openmp" -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX:PATH=$SWEET_LOCAL_SOFTWARE_DST_DIR ../llvm
config_make_default
config_exec make install

# test suits
config_exec make check-clang check-flang

config_success

