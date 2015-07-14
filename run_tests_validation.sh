#! /bin/bash



cd run_tests_validation

for i in *.sh; do
	echo "******************************************************"
	echo "* Executing script $i"
	echo "******************************************************"
	./$i
done
