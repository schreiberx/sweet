#! /bin/bash

for i in run_*.sh; do
	bsub < "$i"
done

