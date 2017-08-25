#! /bin/bash

for i in */run.sh; do
	qsub "$i"
done

