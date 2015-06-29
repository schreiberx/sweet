#! /bin/bash

LIBS="sweet_support"

for i in $LIBS; do
	export LD_LIBRARY_PATH="$HOME/local/$i/lib:$LD_LIBRARY_PATH"
	export LD_LIBRARY_PATH="$HOME/local/$i/lib64:$LD_LIBRARY_PATH"
	export PKG_CONFIG_PATH="$HOME/local/$i/lib/pkgconfig:$PKG_CONFIG_PATH"
	export PATH="$HOME/local/$i/bin:$PATH"
done

