#! /bin/bash

for i in script_*; do
	rm $i/prog_h_*
	rm $i/prog_u_*
	rm $i/prog_v_*
	rm $i/out_*

	rm $i.out
done
