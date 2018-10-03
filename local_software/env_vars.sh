#! /bin/bash

#
# This is a template to set up environment variables correctly
#
# It assumes that you install all your libraries in subdirectories in $SWEETROOT/local_software/local
#


#
# START: Some convenient functions
#
# echo_info [message]:
#	Output in regular colors
#
# echo_success [message]:
#	Output success message in green color
#
# echo_warning [message]:
#	Output warning message in yellow color
#
# echo_error [message]:
#	Output success message in red color
#
# echo_hline [message]:
#	Output a horizontal separator line
#

ECHO_PREFIX="echo -n \"SWEET: \""

# Pretty output
echo_info()( eval ${ECHO_PREFIX}; echo "${@}"; )
echo_success()( eval ${ECHO_PREFIX}; echo -en "\033[0;32m"; echo "${@}"; echo -en "\033[0m"; )
echo_warning()( eval ${ECHO_PREFIX}; echo -en "\033[0;33m"; echo "${@}"; echo -en "\033[0m"; )
echo_error()( eval ${ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; )

# hlines
echo_hline()( echo "SWEET: *************************************************************************"; )
echo_info_hline()( echo_info "*************************************************************************"; )
echo_success_hline()( echo_success "*************************************************************************"; )
echo_warning_hline()( echo_warning "*************************************************************************"; )
echo_error_hline()( echo_error "*************************************************************************"; )

# output error and exit
echo_error_exit(){ echo_error_hline; eval ${ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; echo_error_hline; exit 1; }

#
# END CONV
#


if [ "#$SWEET_ROOT" != "#" ]; then

	echo_warning "Environment variables already loaded (skipping)"
	if [ "`basename -- "$0"`" = "env_vars.sh" ]; then
		return
	fi
else

	echo_hline
	if [ "#$0" != "#-bash" ]; then
		if [ "`basename -- "$0"`" = "env_vars.sh" ]; then
			if [ "`basename -- "$0"`" != "bash" ]; then
				if [ "`basename -- "$0"`" != "modules_env_yellowstone.inc" ]; then
					echo_error ""
					echo_error ">>> $0"
					echo_error "THIS SCRIPT MAY NOT BE EXECUTED, BUT INCLUDED IN THE ENVIRONMENT VARIABLES!"
					echo_error "Use e.g. "
					echo_error ""
					echo_error "   $ source ./env_vars.sh"
					echo_error ""
					echo_error "to setup the environment variables correctly"
					echo_error ""
					exit 1
				fi
			fi
		fi
	fi


	if [ "`basename -- "$SHELL"`" != "bash" ]; then
		echo_error ""
		echo_error "These scripts are only compatible to the bash shell"
		echo_error ""
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
				echo_error "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo_error "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo_error "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo_error ""
				echo_error "Ubuntu system detected and packages missing, please use"
				echo_error ""
				echo_error "    sudo apt-get install $SWEET_SYSTEM_PACKAGES"
				echo_error ""
				echo_error "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo_error "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				echo_error "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
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

	echo_info " + Setting up platform independent environment variables..."

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

	echo_info " + Loading platform specific environment variables..."


	#
	# Include cluster-specific scripts
	#
	if [ "#$SWEET_PLATFORM" != "#" ]; then
		echo_info "                   Mode: SWEET_PLATFORM (=${SWEET_PLATFORM})"
		PLATFORM_ENV_VARS=$(eval echo "${SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM}/env_vars.sh")

		if [ ! -e "$PLATFORM_ENV_VARS" ]; then
			echo_error "File '${PLATFORM_ENV_VARS}' not found!"
			return 1
		fi

		source "$PLATFORM_ENV_VARS" || return 1

	elif [ "#${1}" = "#" ] || [ "#${1}" = "#FORCE" ]; then
		echo_info "                   Mode: Autodetect"

		#
		# Use python code to detect hardware.
		#
		# Note, that this is not required on cluster nodes,
		# since SWEET_PLATFORM will be used to avoid this autodetection.
		#
		function load_this_platform {

			if [ ! -e "SWEETPlatformAutodetect.py" ]; then
				echo_warning "Warning: 'SWEETPlatformAutodetect.py' not found in ${PWD}"
				return 1
			fi
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
				export SWEET_PLATFORM="$p"
				source "./env_vars.sh" || return 1
				break
			fi
		done
	else
		echo_info "                   Mode: Platform environment file provided"

		# Load platform environment variables if specified
		echo_info "     Platform env. file: \$1='${1}'"

		# Set SWEET Platform override
		export SWEET_PLATFORM="$1"
		source "$1" || return 1
	fi

	if [ -z "${SWEET_PLATFORM}" ]; then
		echo_error "INTERNNAL ERROR: No platform detected!!!"
		return
	fi

	export SWEET_PLATFORM_DIR="$(eval echo "${SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM}/")"

	echo_info "         Using platform: '${SWEET_PLATFORM}'"
	echo_info "     Platform directory: '${SWEET_PLATFORM_DIR}'"

	# Back to local software
	cd "$SCRIPTDIR"

	export PS1="[SWEET.$SWEET_PLATFORM] $PS1"


	cd "$BACKDIR"

	echo_success " Environment setup successful (I hope so...)"
	echo_hline

fi
