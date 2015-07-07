#! /bin/bash

PACKAGE_DIRS_AND_FILES="benchmarks data doc python_mods src"
PACKAGE_DIRS_AND_FILES+="
env_vars.sh
INSTALL
Makefile
README
SConstruct
"
 

PACKAGE_NAME="SWEET_`date +%Y_%m_%d`"
PACKAGE_DIR="$PACKAGE_NAME"
PACKAGE_TARBALL="$PACKAGE_NAME.tar.bz2"

echo "Creating package $PACKAGE_NAME"
rm -f -r "$PACKAGE_DIR"
mkdir "$PACKAGE_DIR"

echo " + copying files"
for file in $PACKAGE_DIRS_AND_FILES; do
	cp -r "../$file" "$@" "$PACKAGE_DIR"
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
