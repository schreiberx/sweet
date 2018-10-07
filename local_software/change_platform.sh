
if [ -z "$TMP_SWEET_ROOT" ]; then
	echo "********************************************************************************"
	echo "* ERROR: SWEET environment variables missing!"
	echo "********************************************************************************"
	return 2>/dev/null
	exit 1
fi

if [ -z "$1" ]; then
	echo "Usage:"
	echo "	$0 [platform name]"
	echo
fi

SWEET_PLATFORM="$1"

PLATFORM_ENV_VARS=$(eval echo "${TMP_SWEET_ROOT}/platforms/"??"_${SWEET_PLATFORM}/env_vars.sh")

if [ ! -e "$PLATFORM_ENV_VARS" ]; then
	echo_error "File '${PLATFORM_ENV_VARS}' not found!"
	return 1
fi

source "$PLATFORM_ENV_VARS" || return 1
