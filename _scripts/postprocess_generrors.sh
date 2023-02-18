#! /bin/bash


THISDIR=`pwd`

cd "../../../"

source ./local_software/env_vars.sh || exit 1

SWEETDIR=`pwd`

cd "$THISDIR"

#Verify arguments
if [ $# -eq 0  ] ; then
    echo "Please enter arguments:"
    echo " Reference folder (without /)"
    echo " output file name (optional, e.g. errors.txt)"
    echo " Time (optional, e.g. 00000086400)"
    echo " Variable (optional, e.g. prog_h_pert)"
    echo " Difusion (optional, e.g. 100000000 )"
    exit 0
fi  

ref=$1
if [ -z "$ref" ]; then 
    echo "Please enter arguments:"
    echo " Reference folder (without /)"
    echo " Time (optional, e.g. 00000086400.00000000)"
    echo " variable (optional, e.g. prog_h_pert)"
    echo " output file name (optional, e.g. errors.txt)"
    exit 0
fi
echo $ref

time=$2
if [ -z "$time" ]; then
    time="00000086400.00000000" #1 day
fi
echo $time

var=$3
if [ -z "$var" ]; then
    var="prog_h_pert" # h
fi
echo $var

dif=$4
if [ -z "$dif" ]; then
    dif="100000000" # h
fi
echo $dif

out=$5
if [ -z "$out" ]; then
    out="errors.txt" # output filename
    timestr=${time%.00000000}
    out="Errors_""$var""_t""$timestr""_dif""$dif"".txt"
fi
echo $out


echo ""
echo "Variable Method1 Method1Paper Time dt Variable Method2 Method2Paper Time dt L1Error RMSError MaxError" > $out
echo "Variable Method1 Method1Paper Time dt Variable Method2 Method2Paper Time dt L1Error RMSError MaxError" 


file="output_""$var""_t""$time"".csv"
reffile="$ref""/""$file"

DIRS="script_*u""$dif""_*"

for i in $DIRS; do
	test -d "$i" || continue
	#echo "$i"
	#cd "$i"
	datafile="$i""/""$file"
	#echo $datafile
	#echo "$datafile"
	output="$(python3 ./pp_compute_max_and_rms_errors_interpol.py "$reffile" "$datafile")" #>> "$out"
	echo $output
	if [ "${output:0:1}" != ":" ];then
		echo $output >> $out
	fi
	#errors=`tail -1 "$out"`
	#echo "$i" 
	#echo "$errors"
	#./run.sh | tee "../$i.out"
	#test ${PIPESTATUS[0]} -eq 0 || exit 2>&1
	#cd ".."
done
