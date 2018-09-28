

if [ "${HOSTNAME:0:8}" == "cheyenne" ]; then
#       echo "Loading GNU 7.1.0 module on Cheyenne"
#       module load gnu/7.1.0
	echo "Leaving compiler to mpiCC = intel"
	echo "Loading newer version of git"
	module load git
fi

