#
# Helper functions and environment variables for install scripts
#

source ./env_vars.sh "" || exit 1


#
# Setup variables where to download source, to compile and to install
#
PWD=`pwd`
SWEET_LOCAL_SOFTWARE_SRC_DIR="$PWD/local_src"
SWEET_LOCAL_SOFTWARE_DST_DIR="$PWD/local"

#
# Prepare directory structure
#
mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR"
mkdir -p "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

#
# Determine number of cores to involve in compilation process
#
NPROCS="$(nproc --all)"
if [ "$NPROCS" -gt "10" ]; then
	# We limit the number of parallel build processes
	# This is important on architectures such as Cheyenne where this
	# results in compilation errors due to a lack of resources
	NPROCS=10
fi
MAKE_DEFAULT_OPTS=" -j ${NPROCS}"



function download() {
	echo_info "Downloading from '${1}'"

	#curl -C - "$1" -o "$2" || exit 1
	wget --continue --progress=bar "$1" -O "$2" || exit 1
}


function config_package()
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
	echo_info_hline
	echo_info "Package: ${PKG_NAME}"
	echo_info_hline

	if [ "$1" == "FORCE"  ]; then
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
			exit 0
		fi
	fi

	if [ "${PKG_URL_SRC:0:4}" != "http" ]; then
		PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/${PKG_URL_SRC}"
	fi

	PKG_FILENAME="$(basename ${PKG_URL_SRC})"


	echo_info "Changing to directory '${SWEET_LOCAL_SOFTWARE_SRC_DIR}'"
	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	# Download file
	echo_info "Downloading '${PKG_URL_SRC}'"
	download "$PKG_URL_SRC" "$PKG_FILENAME"

	# Determine how to extract compressed file
	EXT="${PKG_FILENAME##*.}"
	EXTRACT_PROG=""
	if [ "#$EXT" = "#gz" ]; then
		EXTRACT_PROG="tar xzf ${PKG_FILENAME}"
	elif [ "#$EXT" = "#tgz" ]; then
		EXTRACT_PROG="tar xzf ${PKG_FILENAME}"
	elif [ "#$EXT" = "#xz" ]; then
		EXTRACT_PROG="tar xf ${PKG_FILENAME}"
	elif [ "#$EXT" = "#bz2" ]; then
		EXTRACT_PROG="tar xjf ${PKG_FILENAME}"
	else
		echo_error_exit "Unknown extension '${EXT}'"
	fi

	echo_info "Extracting '${PKG_FILENAME}'"
	$EXTRACT_PROG || exit 1

	if [ -z "${PKG_SRC_SUBDIR}" ]; then
		# Automatically detect basename
		PKG_SRC_SUBDIR="${PKG_FILENAME}"
		PKG_SRC_SUBDIR="${PKG_SRC_SUBDIR%%.tar.*}"
		PKG_SRC_SUBDIR="${PKG_SRC_SUBDIR%%.tgz}"
		echo_info "Assuming source foldername to be '${PKG_SRC_SUBDIR}'"
	fi

	if [ ! -e "${PKG_SRC_SUBDIR}" ]; then
		echo_error_exit "Source folder '${PKG_SRC_SUBDIR}' not found"
	fi

	echo_info "Changing to package directory '${PKG_SRC_SUBDIR}'"
	cd "$PKG_SRC_SUBDIR"

	echo_info "Suggesting to use '${NPROCS}' parallel build processes"
}


function config_configure()
{
	M="./configure --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR $@"
	echo_info "Executing '${M}'"
	$M || echo_error_exit "FAILED '${M}'"
}
function config_make_default()	# make
{
	M="make $MAKE_DEFAULT_OPTS $@"
	echo_info "Executing '${M}'"
	$M || echo_error_exit "FAILED '${M}'"
}
function config_make_clean()	# make clean
{
	M="make $MAKE_DEFAULT_OPTS $@ clean"
	echo_info "Executing '${M}'"
	$M || echo_error_exit "FAILED '${M}'"
}
function config_make_install()	# make install
{
	M="make $MAKE_DEFAULT_OPTS $@ install"
	echo_info "Executing '${M}'"
	$M || echo_error_exit "FAILED '${M}'"
}

# Combined functions
function config_make_default_install()			# make; make install
{
	config_make_default
	config_make_install
}
function config_configure_make_default_install()	# ./configure...; make ...; make install
{
	config_configure
	config_make_default
	config_make_default_install
}



function config_success()
{
	echo_success_hline
	echo_success "Package '${PKG_NAME}' successfully installed"
	echo_success_hline
}
