#! /usr/bin/env bash

lscpu -e=CPU,CORE | tail -n +2 | uniq -s 4 -u | sed "s/ .*//"

