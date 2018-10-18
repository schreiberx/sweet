

CHANGE_PLATFORM_PWD=$PWD

if [ -z "$SWEET_ROOT" ]; then
	echo "********************************************************************************"
	echo "* ERROR: SWEET environment variables missing!"
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

	PLATFORM_ENV_VARS=$(eval echo "${SWEET_ROOT}/platforms/"??"_$1/env_vars.sh")

	if [ ! -e "$PLATFORM_ENV_VARS" ]; then
		echo_error "File '${PLATFORM_ENV_VARS}' not found!"
		return 2>/dev/null
		exit 1
	fi

	source "$PLATFORM_ENV_VARS" || return 1

	export SWEET_PLATFORM_ID="$1"

	export SWEET_PLATFORM_DIR="$(eval echo "${SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM_ID}/")"

elif [[ "#$SWEET_PLATFORM_ID" != "#" && "#$1" != "#auto" ]]; then
	#
	# Environment variable "$SWEET_PLATFORM_ID"
	#

	echo_info_hline
	echo_info "                   Mode: manual: SWEET_PLATFORM_ID (=${SWEET_PLATFORM_ID})"
	echo_info_hline

	export SWEET_PLATFORM_DIR="$(eval echo "${SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM_ID}/")"

else

	echo_info_hline
	echo_info "                   Mode: Autodetect"
	echo_info_hline

	#
	# Use python code to detect hardware.
	#
	# Note, that this is not required on cluster nodes,
	# since SWEET_PLATFORM_ID will be used to avoid this autodetection.
	#
	function load_this_platform {

		if [ ! -e "SWEETPlatformAutodetect.py" ]; then
			echo_warning "Warning: 'SWEETPlatformAutodetect.py' not found in ${PWD}"
			return 1
		fi
		# Automagic detection here if called from terminal
		echo -en "import sys\nimport SWEETPlatformAutodetect\nsys.exit(0 if SWEETPlatformAutodetect.autodetect() else 1)" | /usr/bin/env python && return 0

		# Platform not detected
		return 1
	}

	#
	# Try with autodetection of platforms
	#
	for i in $SWEET_ROOT/platforms/??_*; do
		echo "Testing for platform directory $i"
		cd "$i"
		load_this_platform
		if [ $? -eq 0 ]; then
			p=$(basename "$i")
			p=${p/??_/}
			export SWEET_PLATFORM_ID="$p"
			source "$i/env_vars.sh" || return 1
			break
		fi
	done
fi


cd "$CHANGE_PLATFORM_PWD"

if [ -z "${SWEET_PLATFORM_ID}" ]; then
	echo_error_hline
	echo_error "INTERNNAL ERROR: No platform detected!!!"
	echo_error_hline
	return
fi

export SWEET_PLATFORM_DIR="$(eval echo "${SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM_ID}/")"

echo_info "         Using platform: '${SWEET_PLATFORM_ID}'"
echo_info "     Platform directory: '${SWEET_PLATFORM_DIR}'"


echo_success_hline
echo_success "Platform changed to '${SWEET_PLATFORM_ID}'"
echo_success_hline
