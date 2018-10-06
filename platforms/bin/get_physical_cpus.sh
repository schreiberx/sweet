#! /usr/bin/env bash

lscpu -e=CPU,CORE | tail -n +2 | sort -k2,2n -k1,1n | uniq -s 3 | sed "s/ .*//"

