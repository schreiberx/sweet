#! /bin/bash


BASENAME="`pwd`"
for i in script_*out; do 
	TAG="Simulation time (seconds):"
	echo "$i" $(grep "$TAG" "$i" | sed "s/$TAG//")
done
