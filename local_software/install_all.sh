#! /bin/bash

source ./config_install.sh "" || exit 1
source ./env_vars.sh "" || exit 1


if [ ! -e "$SWEET_PLATFORM_DIR" ]; then
	echo_error ""
	echo_error "SWEET Platform directory '${SWEET_PLATFORM_DIR}' not found"
	echo_error ""
	exit 1
fi

LOCAL_SOFTWARE_PLATFORM="${SWEET_PLATFORM_DIR}/local_software_default.sh"
if [ ! -e "$LOCAL_SOFTWARE_PLATFORM" ]; then
	echo_error ""
	echo_error "ERROR: File '${LOCAL_SOFTWARE_PLATFORM}' not found"
	echo_error ""
	exit
fi


# Packages to install
PKGS=()
source "$LOCAL_SOFTWARE_PLATFORM"


echo_info "Installing packages: ${PKGS[@]}"

for D in "${PKGS[@]}"; do
	./$D || exit 1
done

echo_hline
echo_success "All software installed"
echo_hline

