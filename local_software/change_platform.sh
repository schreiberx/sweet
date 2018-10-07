
if [ -z "$TMP_SWEET_ROOT" ]; then
	echo "********************************************************************************"
	echo "* ERROR: SWEET environment variables missing!"
	echo "********************************************************************************"
	return 2>/dev/null
	exit 1
fi

if [ -z "$1" ]; then
	echo_info "Usage:"
	echo_info "	$0 [platform name]"
	echo_info
fi

PLATFORM_ENV_VARS=$(eval echo "${TMP_SWEET_ROOT}/platforms/"??"_$1/env_vars.sh")

if [ ! -e "$PLATFORM_ENV_VARS" ]; then
	echo_error "File '${PLATFORM_ENV_VARS}' not found!"
	return 2>/dev/null
	exit 1
fi

source "$PLATFORM_ENV_VARS" || return 1

export SWEET_PLATFORM="$1"

export SWEET_PLATFORM_DIR="$(eval echo "${TMP_SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM}/")"

echo_success_hline
echo_success "Platform changed to '${SWEET_PLATFORM}'"
echo_success_hline
