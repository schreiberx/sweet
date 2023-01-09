#
# Helper functions and environment variables for install scripts
#


if [ -z "$MULE_SOFTWARE_ROOT" ]; then
	echo "****************************************************************************************************"
	echo "* ERROR: SWEET environment variables missing!"
	echo "****************************************************************************************************"
	exit 1
fi


#
# Setup variables where to download source, to compile and to install
#
PWD_=`pwd`

SWEET_LOCAL_SOFTWARE_SRC_DIR="$PWD_/local_src"
SWEET_LOCAL_SOFTWARE_DST_DIR="$PWD_/local"


#
# Prepare directory structure
#
mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR"
mkdir -p "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

#
# Determine number of cores to involve in compilation process
#
if [ "`uname`" == "Darwin" ]; then
	MAKE_DEFAULT_OPTS="-j"

else
	NPROCS="$(nproc)"
	if false; then
		if [ "$NPROCS" -gt "10" ]; then
			# We limit the number of parallel build processes
			# This is important on architectures such as Cheyenne where this
			# results in compilation errors due to a lack of resources
			NPROCS=10
		fi
	fi

	MAKE_DEFAULT_OPTS=" -j ${NPROCS}"
fi


#
# Travis build server requires output every 10 minutes, otherwise it's assumed to be dead and killed
#
# We use this alive thread to output some dots
#
function config_alive_thread()
{
	PREV_OUTPUT_LINES=$(wc -l < "${PKG_CONFIG_STD_OUTPUT}")
	while true; do
		OUTPUT_LINES=$(wc -l < "${PKG_CONFIG_STD_OUTPUT}")
		C=$((OUTPUT_LINES-PREV_OUTPUT_LINES))
		while [[ $C -gt 0 ]]; do
			echo -n "."
			C=$(($C-1))
		done
		PREV_OUTPUT_LINES=$OUTPUT_LINES
		sleep 1
	done
}

function config_alive_thread_start()
{
	config_alive_thread &
	CONFIG_ALIVE_THREAD_PID=$!

	echo -n "MULE: Status ["
}

function config_alive_thread_stop()
{
	kill -0 $CONFIG_ALIVE_THREAD_PID
	if [[ $? -eq 0 ]]; then
		# Disown thread to avoid any "Terminated" message
		disown $CONFIG_ALIVE_THREAD_PID
		kill $CONFIG_ALIVE_THREAD_PID >/dev/null 2>&1
		echo "] done"
	fi
}


function config_setup()
{
	if [[ -z "$PKG_NAME" ]]; then
		echo_error_exit "config_setup() PKG_NAME is not set"
	fi

	echo_info_hline
	echo_info "Package: ${PKG_NAME}"
	echo_info_hline

	if [[ ! -z "DEBUG" ]]; then
		export PKG_CONFIG_STD_OUTPUT="${SWEET_LOCAL_SOFTWARE_SRC_DIR}/${PKG_NAME}_config.out"

		echo_info "Writing installer output to '${PKG_CONFIG_STD_OUTPUT}'"
		echo -n "" > "$PKG_CONFIG_STD_OUTPUT"
	else
		export PKG_CONFIG_STD_OUTPUT="/dev/stdout"
	fi
}



function config_error_exit()
{
	# Kill alive thread if it exists
	config_alive_thread_stop

	echo_info "ERROR occured, please check log file in ./local_src"

	#
	# print output message
	# print log file
	# exit
	#
	echo_error_hline
	echo_error "ERROR: $@"
	echo_error_hline
	echo_error "An error occurred during the build"
	echo_error "Showing last 100 lines from build output file"
	echo_error "$PKG_CONFIG_STD_OUTPUT"

	tail -n 100 "$PKG_CONFIG_STD_OUTPUT"
	#cat "$PKG_CONFIG_STD_OUTPUT"
	exit 1
}



function config_extract_fun()
{
	EXT="${0##*.}"
	EXTRACT_PROG=""
	if [ "#$EXT" = "#gz" ]; then
		TAR_CMD="zvf"
	elif [ "#$EXT" = "#tgz" ]; then
		TAR_CMD="zvf"
	elif [ "#$EXT" = "#xz" ]; then
		TAR_CMD="fv"
	elif [ "#$EXT" = "#bz2" ]; then
		TAR_CMD="jvf"
	elif [ "#$EXT" = "#sh" ]; then
		echo "Shell script detected, skipping extraction"
		return
	else
		config_error_exit "Unknown extension '${EXT}'"
	fi

	EXTRACT_PROG="tar x${TAR_CMD} ${PKG_FILENAME}"
	LIST_TAR="tar t${TAR_CMD} ${PKG_FILENAME}"
	echo $LIST_TAR

	echo_info "Extracting '${EXTRACT_PROG}'"
	$EXTRACT_PROG || config_error_exit 'Failed to extract archive'
}

function config_extract()
{
	echo "Exctracting $@"
	config_alive_thread_start
	config_extract_fun $@ >> $PKG_CONFIG_STD_OUTPUT 2>&1
	config_alive_thread_stop
}



function config_download_fun()
{
	PKG_FILENAME="$(basename ${1})"

	if [ "${MULE_PLATFORM_ID:0:9}" == "supermuc_" ]; then
		echo_warning "Detected SUPERMUC platform, skipping download!" 1>&2

		if [ ! -f "${PKG_FILENAME}" ]; then
			echo_error "Skipped download on SUPERMUC, but source file '${PKG_FILENAME}' not found."
			echo_error "Please copy it manually to the local_src folder!"
			exit 1
		fi

		return
	fi

        if type wget >/dev/null 2>&1; then
		echo_info "Using 'wget'" 1>&2
		wget --continue --progress=bar "$1" -O "$PKG_FILENAME" || config_error_exit "Download failed! Did you install the certificates via ./install_cacerts.sh"
	else
		echo_info "Using 'curl'" 1>&2
		# Do not continue download since some curl versions seem to be buggy :-(
		curl -o "$PKG_FILENAME" "$1" || config_error_exit "Download failed! Did you install the certificates via ./install_cacerts.sh"
	fi
}

function config_download()
{
	echo_info "Downloading $@"
	config_download_fun $@ >> $PKG_CONFIG_STD_OUTPUT #2>&1
}



function config_download_and_extract()
{
	config_download $@
	config_extract $@
}



function config_package_test_existing_dir_fun()
{
	#
	# PKG_NAME: string
	#	Name of package. This is just for convenience.
	#
	# PKG_INSTALLED_FILE: string
	#	Path to file which would exist if this package is installed
	#
	# PKG_URL_SRC: string
	#	URL to source.
	#	If the URL doesn't start with http, the prefix
	#		https://www.martin-schreiber.info/pub/sweet_local_software/
	#	is used for sake of convenience
	#
	# PKG_SRC_SUBDIR: string
	#	Name of folder in which the source code can be found after extracting the package
	#
	#echo ${PKG_INSTALLED_FILE}
	if [ "$1" = "FORCE"  ]; then
		echo_warning "FORCE detected => reinstallation of package"
	else
		if [ -e "${PKG_INSTALLED_FILE}" ]; then
			echo_warning_hline
			echo_warning "The package '${PKG_NAME}' is already installed"
			echo_warning "Use"
			echo_warning "	${0} FORCE"
			echo_warning "to reinstall package or"
			echo_warning "	${0} FORCE CLEAN"
			echo_warning "to remove existing source files as well."
			echo_warning_hline
			exit 0	# Exit with success code since software was already successfully installed
		fi
	fi
}

function config_package_test_existing_dir()
{
	config_package_test_existing_dir_fun $@	# We really want to get the message that this package was already installed
						# >> "$PKG_CONFIG_STD_OUTPUT" 2>&1
}

function config_package_download()
{
	if [ "${PKG_URL_SRC:0:4}" != "http" ]; then
		export PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/${PKG_URL_SRC}"
	fi

	echo_info "Changing to directory '${SWEET_LOCAL_SOFTWARE_SRC_DIR}'"
	config_exec cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	config_download  "$PKG_URL_SRC"
}


function config_package_extract()
{
	if [[ -z "$PKG_FILENAME" ]]; then
		echo_error "Package filename not set, skipping extraction but generating directory"

		echo_info "Changing to directory '${SWEET_LOCAL_SOFTWARE_SRC_DIR}'"
		cd "${SWEET_LOCAL_SOFTWARE_SRC_DIR}"

		PKG_SRC_SUBDIR="$(echo -n "${PKG_NAME}" | tr '[:upper:]' '[:lower:]')"
		mkdir -p "${PKG_SRC_SUBDIR}"
		cd "${PKG_SRC_SUBDIR}"
		return
	fi

	# Determine how to extract compressed file
	# !!!
	# Be verbose to determine directory!!!
	# !!!
	EXT="${PKG_FILENAME##*.}"
	EXTRACT_PROG=""
	if [ "#$EXT" = "#gz" ]; then
		TAR_CMD="zf"
	elif [ "#$EXT" = "#tgz" ]; then
		TAR_CMD="zf"
	elif [ "#$EXT" = "#xz" ]; then
		TAR_CMD="f"
	elif [ "#$EXT" = "#bz2" ]; then
		TAR_CMD="jf"
	elif [ "#$EXT" = "#sh" ]; then
		echo "Shell script detected, skipping extraction"
		return
	else
		config_error_exit "Unknown extension '${EXT}'"
	fi

	EXTRACT_PROG="tar x${TAR_CMD} ${PKG_FILENAME}"
	LIST_TAR="tar t${TAR_CMD} ${PKG_FILENAME}"

	if [ -z "${PKG_SRC_SUBDIR}" ]; then
		echo_info "Detecting target directory"
		LIST_OUTPUT="$($LIST_TAR || config_error_exit 'Failed to determine content of archive')"
		PKG_SRC_SUBDIR=$(echo "$LIST_OUTPUT" | head -n 1 | sed "s/\/.*//")
		echo_info "Detected source foldername to be '${PKG_SRC_SUBDIR}'"
	fi

	if [ "$2" == "CLEAN"  ]; then
		echo_warning "CLEAN detected => removing source folder '${PKG_SRC_SUBDIR}'"
		if [ ! -e "${PKG_SRC_SUBDIR}" ]; then
			echo_warning "Directory '${PKG_SRC_SUBDIR}' does not exist"
		else
			rm -rf "${PKG_SRC_SUBDIR}"
		fi
	fi

	echo_info "Extracting '${EXTRACT_PROG}'"
	config_alive_thread_start
	$EXTRACT_PROG || config_error_exit 'Failed to extract archive'
	config_alive_thread_stop

	if [ ! -e "${PKG_SRC_SUBDIR}" ]; then
		config_error_exit "Source folder '${PKG_SRC_SUBDIR}' not found"
	fi

	echo_info "Changing to package directory '${PKG_SRC_SUBDIR}'"
	cd "$PKG_SRC_SUBDIR"
}



function config_configure()
{
	M="./configure --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR $@"
	echo_info "Executing '${M}'"

	config_alive_thread_start
	$M >> "$PKG_CONFIG_STD_OUTPUT" 2>&1 || config_error_exit "FAILED '${M}'"
	config_alive_thread_stop
}
function config_make_default()	# make
{
	M="make $MAKE_DEFAULT_OPTS $@"
	echo_info "Executing '${M}'"

	config_alive_thread_start
	$M >> "$PKG_CONFIG_STD_OUTPUT" 2>&1 || config_error_exit "FAILED '${M}'"
	config_alive_thread_stop
}
function config_make_clean()	# make clean
{
	M="make $MAKE_DEFAULT_OPTS $@ clean"
	echo_info "Executing '${M}'"

	config_alive_thread_start
	$M >> "$PKG_CONFIG_STD_OUTPUT" 2>&1 || config_error_exit "FAILED '${M}'"
	config_alive_thread_stop
}
function config_make_install()	# make install
{
	M="make $MAKE_DEFAULT_OPTS $@ install"
	echo_info "Executing '${M}'"

	config_alive_thread_start
	$M >> "$PKG_CONFIG_STD_OUTPUT" 2>&1 || config_error_exit "FAILED '${M}'"
	config_alive_thread_stop
}
function config_exec()
{
	M="$@"
	echo_info "Executing '${M}'"

	config_alive_thread_start
	"${@}" >> "$PKG_CONFIG_STD_OUTPUT" 2>&1 || config_error_exit "FAILED '${M}'"
	config_alive_thread_stop
}


function config_package()
{
	echo_info "Testing for existing package installation..."
	config_package_test_existing_dir $@

	echo_info "Downloading package..."
	config_package_download $@

	echo_info "Extracting package..."
	config_package_extract $@

	echo_info "Suggesting to use '${NPROCS}' parallel build processes"
}


# Combined functions
function config_configure_make_default()	# ./configure...; make ...;
{
	config_configure ${@}
	config_make_default
}
function config_make_default_install()			# make; make install
{
	config_make_default ${@}
	config_make_install
}
function config_configure_make_default_install()	# ./configure...; make ...; make install
{
	config_configure ${@}
	config_make_default
	config_make_install
}



function config_success()
{
	echo_success_hline
	echo_success "Package '${PKG_NAME}' successfully installed"
	echo_success_hline
}
