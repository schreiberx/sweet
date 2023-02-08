

CHANGE_PLATFORM_PWD=$PWD

if [ -z "$MULE_ROOT" ]; then
	echo "********************************************************************************"
	echo "* ERROR: MULE_ROOT environment variables missing!"
	echo "********************************************************************************"
	return 2>/dev/null
	exit 1
fi

#if [ -z "$1" ]; then
#	echo_info "Usage:"
#	echo_info "	${BASH_SOURCE[0]} [platform name]"
#	echo_info
#	return
#fi


if [[ "#$1" != "#" && "#$1" != "#auto" ]]; then

	#
	# Parameter of this script
	#
	echo_info_hline
	echo_info "                   Mode: manual: $0 $1 (=${1})"
	echo_info_hline

	PLATFORM_ENV_VARS=$(eval echo "${MULE_ROOT}/platforms/"??"_$1/env_vars.sh")

	if [ ! -e "$PLATFORM_ENV_VARS" ]; then
		PLATFORM_ENV_VARS=$(eval echo "${MULE_ROOT}/platforms/$1/env_vars.sh")

		if [ ! -e "$PLATFORM_ENV_VARS" ]; then
			echo_error "File '${PLATFORM_ENV_VARS}' not found!"
			return 2>/dev/null
			exit 1
		fi

		export MULE_PLATFORM_ID="${1:3}"
	else
		export MULE_PLATFORM_ID="${1}"
	fi

	source "$PLATFORM_ENV_VARS"

	export MULE_PLATFORM_DIR="$(eval echo "${MULE_ROOT}/platforms/"??"_${MULE_PLATFORM_ID}/")"

elif [[ "#$MULE_PLATFORM_ID" != "#" && "#$1" != "#auto" ]]; then
	#
	# Environment variable "$MULE_PLATFORM_ID"
	#

	echo_info_hline
	echo_info "                   Mode: manual: MULE_PLATFORM_ID (=${MULE_PLATFORM_ID})"
	echo_info_hline

	export MULE_PLATFORM_DIR="$(eval echo "${MULE_ROOT}/platforms/"??"_${MULE_PLATFORM_ID}/")"

else

	echo_info_hline
	echo_info "                   Mode: Autodetect"
	echo_info_hline

	#
	# Use python code to detect hardware.
	#
	# Note, that this is not required on cluster nodes,
	# since MULE_PLATFORM_ID will be used to avoid this autodetection.
	#
	function load_this_platform {

		if [ ! -e "JobPlatformAutodetect.py" ]; then
			echo_warning "Warning: 'JobPlatformAutodetect.py' not found in ${PWD}"
			return 1
		fi
		# Automagic detection here if called from terminal
		echo -en "import sys\nimport JobPlatformAutodetect\nsys.exit(0 if JobPlatformAutodetect.autodetect() else 1)" | /usr/bin/env python3 && return 0

		# Platform not detected
		return 1
	}

	#
	# Try with autodetection of platforms
	#
	#echo "Testing debug output " $MULE_ROOT/platforms/??_*
	for i in $MULE_ROOT/platforms/??_*; do
		# Test whether it's really a proper platform directory and no leftover (including __pycache__)
		if [[ ! -e "$i/JobPlatform.py" ]]; then
			echo "$i/JobPlatform.py not found, but platform directory exists"
			continue
		fi

		cd "$i"
		load_this_platform
		if [ $? -eq 0 ]; then
			p=$(basename "$i")
			p=${p/??_/}
			export MULE_PLATFORM_ID="$p"
			source "$i/env_vars.sh" || return 1
			break
		fi
	done
fi


cd "$CHANGE_PLATFORM_PWD"

if [ -z "${MULE_PLATFORM_ID}" ]; then
	echo_error_hline
	echo_error "INTERNNAL ERROR: No platform detected!!!"
	echo_error_hline
	exit 1
fi

export MULE_PLATFORM_DIR="$(eval echo "${MULE_ROOT}/platforms/"??"_${MULE_PLATFORM_ID}/")"
export MULE_LOCAL_PLATFORM_DIR="$(eval echo "${MULE_LOCAL_ROOT}/platforms/"??"_${MULE_PLATFORM_ID}/")"

echo_info "               Using platform: '${MULE_PLATFORM_ID}'"
echo_info "      Platform mule directory: '${MULE_PLATFORM_DIR}'"
echo_info "Platform mule_local directory: '${MULE_LOCAL_PLATFORM_DIR}'"


echo_success_hline
echo_success "Platform changed to '${MULE_PLATFORM_ID}'"
echo_success_hline
