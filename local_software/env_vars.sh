#! /bin/bash

#
# This is a template to set up environment variables correctly
#
# It assumes that you install all your libraries in subdirectories in $SWEETROOT/local_software/local
#


if [ "#$SWEET_ROOT" != "#" ]; then
	echo "SWEET environment variables already loaded (skipping)"
	if [ "`basename -- "$0"`" == "env_vars.sh" ]; then
		return
	fi
else

	if [ "#$0" != "#-bash" ]; then
		if [ "`basename -- "$0"`" == "env_vars.sh" ]; then
			if [ "`basename -- "$0"`" != "bash" ]; then
				if [ "`basename -- "$0"`" != "modules_env_yellowstone.inc" ]; then
					echo "ERROR|"
					echo "ERROR| >>> $0"
					echo "ERROR| THIS SCRIPT MAY NOT BE EXECUTED, BUT INCLUDED IN THE ENVIRONMENT VARIABLES!"
					echo "ERROR| Use e.g. "
					echo "     |"
					echo "     |    $ source ./env_vars.sh"
					echo "     |"
					echo "ERROR| to setup the environment variables correctly"
					echo "ERROR|"
					exit
				fi
			fi
		fi
	fi


	if [ "`basename -- "$SHELL"`" != "bash" ]; then
		echo "ERROR|"
		echo "ERROR| These scripts are only compatible to the bash shell"
		echo "ERROR|"
		return
	fi


	#
	# Detect Ubuntu system
	#
	grep "Ubuntu" /etc/issue > /dev/null
	if [ "x$?" = "x0" ]; then
		export SWEET_SYSTEM_PACKAGES="libxft-dev libssl-dev"

		#
		# Check for extra package if X-server is locally activated
		#
		if [ "x$DISPLAY" = "x:0" ]; then
			export SWEET_SYSTEM_PACKAGES="$SWEET_SYSTEM_PACKAGES libgl1-mesa-dev"
		fi

		for i in $SWEET_SYSTEM_PACKAGES; do
			dpkg -s "$i" >/dev/null 2>&1
			if [ "x$?" != "x0" ]; then
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo ""
				echo "Ubuntu system detected and packages missing, please use"
				echo ""
				echo "    sudo apt-get install $SWEET_SYSTEM_PACKAGES"
				echo ""
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				return
			fi
		done
	fi

	# Backup current directory
	BACKDIR="$PWD"

	SCRIPTDIR="$(dirname "${BASH_SOURCE[0]}")"
	cd "$SCRIPTDIR"
	# Get full path
	SCRIPTDIR="$PWD"

	# Get SWEET root directory
	cd "../"
	export SWEET_ROOT="$PWD"

	# Back to local software
	cd "$SCRIPTDIR"


	#
	# Include cluster-specific scripts
	#
	if [ "#$SWEET_PLATFORM" != "#" ]; then

		ENV_VARS=$(eval echo "${SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM}/env_vars.sh")
		echo "Loading SWEET_PLATFORM='${SWEET_PLATFORM}' platform environment variables from ${ENV_VARS}"
		source "$ENV_VARS"

	elif [ "#$1" == "#" ]; then
		#
		# Use python code to detect hardware.
		#
		# Note, that this is not required on cluster nodes,
		# since SWEET_PLATFORM will be used to avoid this autodetection.
		#
		function load_this_platform {

			# Automagic detection here if called from terminal
			echo -en "import sys\nimport SWEETPlatformAutodetect\nsys.exit(0 if SWEETPlatformAutodetect.autodetect() else 1)" | /usr/bin/env python3 && return 0

			# Platform not detected
			return 1
		}

		#
		# Try with autodetection of platforms
		#
		for i in $SWEET_ROOT/platforms/??_*; do
			cd "$i"
			load_this_platform
			if [ $? -eq 0 ]; then
				p=$(basename "$i")
				p=${p/??_/}
				echo "SWEET: Detected platform '$p'"
				export SWEET_PLATFORM="$p"
				source "./env_vars.sh"
				break
			fi
		done
	else
		# Load platform environment variables if specified
		source "$1"
	fi

	# Back to local software
	cd "$SCRIPTDIR"


	if [ ! -d "$SCRIPTDIR" ]; then
		echo
		echo "ERROR| Execute this script only from the SWEET root directory"
		echo "     |   $ source local_software/env_vars.sh"
		echo
		cd "$BACKDIR"
		return
	fi

	export PS1="[SWEET] $PS1"

	export PATH="$SCRIPTDIR/local/bin:$PATH"
	export PKG_CONFIG_PATH="$SCRIPTDIR/local/lib/pkgconfig:$PKG_CONFIG_PATH"

	export LD_LIBRARY_PATH="$SCRIPTDIR/local/lib:$LD_LIBRARY_PATH"
	if [ -d "$SCRIPTDIR/local/lib64" ]; then
		export LD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$LD_LIBRARY_PATH"
	fi

	export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib:$LD_LIBRARY_PATH"
	if [ -d "$SCRIPTDIR/local/lib64" ]; then
		export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$LD_LIBRARY_PATH"
	fi

	export PYTHONPATH="$PYTHONPATH:$SCRIPTDIR/local/lib/python3.6/site-packages/"

	# Add SWEET python_mods to pythonpath
	export PYTHONPATH="$PYTHONPATH:$SWEET_ROOT/python_mods"

	echo "SWEET: SWEET environment variables loaded"

	cd "$BACKDIR"


fi
