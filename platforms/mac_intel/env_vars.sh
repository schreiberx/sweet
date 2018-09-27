

if [ "${HOSTNAME:0:15}" == "mac-login-intel" -o "${HOSTNAME:0:7}" == "mac-snb" ]; then
	echo "Loading modules for mac-login-intel"

	echo "Loading GCC/7"
	module unload gcc
	module load gcc/7

	echo "Loading binutils"
	module load binutils/2.25

#       module unload intel
#       module load intel/18.0
fi


