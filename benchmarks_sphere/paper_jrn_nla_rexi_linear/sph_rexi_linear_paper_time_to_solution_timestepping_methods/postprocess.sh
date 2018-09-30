#! /bin/bash


BASENAME="`pwd`"
for i in script_*/output.out; do 
	TAG="Wallclock time (seconds):"
	echo "$i" $(grep "$TAG" "$i" | sed "s/$TAG//")
done
