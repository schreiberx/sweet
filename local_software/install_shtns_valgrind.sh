#! /bin/bash

# AVX512f makes trouble for valgrind
export CFLAGS="$CFLAGS -mno-avx512f"

./install_shtns.sh
