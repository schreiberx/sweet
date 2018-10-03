
#
# Helper functions and environment variables for install scripts
#

PWD=`pwd`
SWEET_LOCAL_SOFTWARE_SRC_DIR="$PWD/local_src"
SWEET_LOCAL_SOFTWARE_DST_DIR="$PWD/local"

mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR"
mkdir -p "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

MAKE_DEFAULT_OPTS=" -j $(nproc)"


function download() {
	echo_info "Downloading from '${1}'"

	#curl -C - "$1" -o "$2" || exit 1
	wget --continue --progress=bar "$1" -O "$2" || exit 1
}


function config_package()
{
	echo_hline
	echo_info "Package: ${PKG_NAME}"
	echo_hline

	if [ "$1" == "FORCE"  ]; then
		echo_warning "FORCE detected => reinstallation of package"
	else
		if [ -e "${PKG_INSTALLED_FILE}" ]; then
			echo_warning_hline
			echo_warning "The package '${PKG_NAME}' is already installed"
			echo_warning "Use"
			echo_warning "	${0} FORCE"
			echo_warning "to reinstall package and"
			echo_warning "	${0} FORCE CLEAN"
			echo_warning "to remove existing source files first."
			echo_warning_hline
			exit 1
		fi
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

	if [ -z "${SRC_SUBDIR}" ]; then
		# Automatically detect basename
		SRC_SUBDIR="${PKG_FILENAME}"
		SRC_SUBDIR="${SRC_SUBDIR%%.tar.*}"
		echo_info "Assuming source foldername to be '${SRC_SUBDIR}'"
	fi

	if [ ! -e "${SRC_SUBDIR}" ]; then
		echo_error_exit "Source folder '${SRC_SUBDIR}' not found"
	fi

	echo_info "Changing to package directory '${SRC_SUBDIR}'"
	cd "$SRC_SUBDIR"
}


function config_make_install()
{
	echo_info "Executing 'make'"
	make $MAKE_DEFAULT_OPTS || echo_error_exit "FAILED make"

	echo_info "Executing 'make install'"
	make $MAKE_DEFAULT_OPTS install || echo_error_exit "FAILED make install"
}


function config_success()
{
	echo_success_hline
	echo_success "Package '${PKG_NAME}' successfully installed"
	echo_success_hline
}
