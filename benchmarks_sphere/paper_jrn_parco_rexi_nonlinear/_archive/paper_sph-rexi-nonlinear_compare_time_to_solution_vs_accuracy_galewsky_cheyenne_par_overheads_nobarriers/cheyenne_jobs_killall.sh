#! /bin/bash

qdel `qselect -u $USER`
#for i in `qstat -u $USER | sed "s/ .*//"`
	
