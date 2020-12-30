#! /usr/bin/env bash

OUTPUT=$(lscpu -p=CPU,CORE || exit 1)
echo "$OUTPUT" | grep "^[^#]" | sort -k2,2n -k1,1n | sed "s/,.*//"

