#! /bin/bash

for i in `bjobs  -u martins | sed "s/ .*//"`; do
	echo "KILLING $i"
	bkill "$i"
done
