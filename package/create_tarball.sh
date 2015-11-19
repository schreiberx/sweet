#! /bin/bash

PACKAGE_DIRS_AND_FILES=" data doc python_mods src run_tests_validation "
#PACKAGE_DIRS_AND_FILES+=" benchmarks "
PACKAGE_DIRS_AND_FILES+="
run_tests_compile_all.sh
run_tests_validation.sh
license.txt
INSTALL
Makefile
README
SConstruct
"

EXCLUDE_FILES="
doc/rexi/results_plots
"
 

PACKAGE_NAME="SWEET_`date +%Y_%m_%d`"
PACKAGE_DIR="$PACKAGE_NAME"
PACKAGE_TARBALL="$PACKAGE_NAME.tar.bz2"

echo "Creating package $PACKAGE_NAME"
rm -f -r "$PACKAGE_DIR"
mkdir "$PACKAGE_DIR"


echo " + copying files"
for file in $PACKAGE_DIRS_AND_FILES; do
	cp -r "../$file" "$PACKAGE_DIR"
done

echo " + handling local software"
#
# local software handler
#
mkdir "$PACKAGE_DIR/local_software"

cp ../local_software/install_* "$PACKAGE_DIR/local_software"
cp ../local_software/clean.sh  "$PACKAGE_DIR/local_software"
cp ../local_software/config.sh "$PACKAGE_DIR/local_software"
cp ../local_software/env_vars.sh "$PACKAGE_DIR/local_software"


echo " + removing excluded files"
for E in $EXCLUDE_FILES; do
	rm -r "$PACKAGE_DIR/$E"
done

#echo " + removing svn information"
## remove svn from package directory
#cd "$PACKAGE_DIR" && { find ./ -name ".svn" | xargs rm -Rf; } && cd ..

echo " + creating tarball $PACKAGE_TARBALL"
rm -f "$PACKAGE_TARBALL"
tar cjf "$PACKAGE_TARBALL" "$PACKAGE_DIR"

echo " + cleaning up"
rm -r "$PACKAGE_DIR"

echo "Done!"
