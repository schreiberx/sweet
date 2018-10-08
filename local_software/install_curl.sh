#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="curl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/curl"
PKG_URL_SRC="curl-7.61.1.tar.xz"

config_package $@ || exit 1

CONFIG_PARAMS=""
CONFIG_PARAMS+=" --enable-http"
CONFIG_PARAMS+=" --enable-file"
CONFIG_PARAMS+=" --enable-proxy"
CONFIG_PARAMS+=" --enable-tftp"
CONFIG_PARAMS+=" --enable-ipv6"

CONFIG_PARAMS+=" --enable-versioned-symbols"

CONFIG_PARAMS+=" --with-nghttp2"
#CONFIG_PARAMS+=" --with-libssh2"

# SSL/TLS config
CONFIG_PARAMS+=" --with-ssl=$SWEET_LOCAL_SOFTWARE_DST_DIR"
CONFIG_PARAMS+=" --with-ca-path=$SWEET_LOCAL_SOFTWARE_DST_DIR/ssl/certs"
CONFIG_PARAMS+=" --with-ca-bundle=$SWEET_LOCAL_SOFTWARE_DST_DIR/ssl/ca-bundle.crt"

config_configure $CONFIG_PARAMS || exit 1

config_make_default || exit 1

if [ "${HOSTNAME:0:10}" != "mpp2-login" ]; then
	# Just save some time :-)
	config_exec "make test" || exit 1
fi
config_make_install || exit 1

echo_info "Installing certificates"
echo_warning "We enforce getting certificates even over an unsecure connection!"
config_exec ./lib/mk-ca-bundle.pl -k

echo_info "Installing bundle"
config_exec cp "ca-bundle.crt" "$SWEET_LOCAL_SOFTWARE_DST_DIR/ssl/ca-bundle.crt" || echo_error_exit "Failed to install ca-bundle.crt"

config_success || exit 1

