#! /bin/bash


THISDIR=`pwd`

#Verify arguments
files=depth*Errors12max.txt
outfile=Errors12max.txt
echo "method/dt Error1 Error2 ErrorMax" > $outfile
for i in $files; do
    echo "$i"
    errors=`cat $i`
    echo $i $errors >> $outfile
done
#cat $outfile
