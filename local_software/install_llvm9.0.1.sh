#! /bin/bash

source ./install_helpers.sh ""

#
# Package based on GIT release package
# https://github.com/llvm/llvm-project/releases/tag/llvmorg-9.0.1
#

PKG_NAME="llvm_9.0.1"
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/llvmorg-9.0.1.tar.gz"

config_setup

config_package $@

# Compile
mkdir -p build
cd build
config_exec cmake -G "Unix Makefiles" ../llvm
config_exec make

# test suits
config_exec make check-clang

config_success

