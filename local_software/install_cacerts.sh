#! /bin/bash

source ./install_helpers.sh "" || exit 1


PKG_NAME="ca-certificates"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/ssl/certs/76579174.0"
#PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/ssl/certs/XRamp_Global_CA_Root.crt"
PKG_URL_SRC="ca-certificates-mozilla_2018_10_05.tar.bz2"

config_package $@

CERTDIR="${SWEET_LOCAL_SOFTWARE_DST_DIR}/ssl/certs"
mkdir -p "${CERTDIR}"

#cd ${CERTDIR} || exit 1

#CAFILE="cacert.pem"
#echo "Splitting $CAFILE"
#csplit -f cert- $CAFILE '/-----BEGIN CERTIFICATE-----/' '{*}'

for i in *; do
	echo "Processing $i"

	config_exec cp "${i}" "${CERTDIR}" || exit 1
	HASH=$(openssl x509 -hash -noout -in "${i}")
	config_exec ln -sf ./$i "${CERTDIR}/${HASH}.0" || exit 1
done
