#! /bin/bash

#
# This is a template to set up environment variables correctly
#
# It assumes that you install all your libraries in subdirectories in $HOME/local
#

if [ "#$0" != "#-bash" ]; then
	if [ "`basename $0`" == "env_vars.sh" ]; then
		if [ "`basename $0`" != "bash" ]; then
			if [ "`basename $0`" != "modules_env_yellowstone.inc" ]; then
				echo "ERROR|"
				echo "ERROR| >>> $0"
				echo "ERROR| THIS SCRIPT MAY NOT BE EXECUTED, BUT INCLUDED IN THE ENVIRONMENT VARIABLES!"
				echo "ERROR| Use e.g. "
				echo "     |"
				echo "     |    $ source ./env_vars.sh"
				echo "     |"
				echo "ERROR| to setup the environment variables correctly"
				echo "ERROR|"
				return
			fi
		fi
	fi
fi

if [ "${HOSTNAME:0:8}" == "cheyenne" ]; then
#	echo "Loading GNU 7.1.0 module on Cheyenne"
#	module load gnu/7.1.0
	echo "Leaving compiler to mpiCC = intel"
	echo "Loading newer version of git"
	module load git
fi


#
# Detect Ubuntu system
#
grep "Ubuntu" /etc/issue > /dev/null
if [ "x$?" = "x0" ]; then
	export SWEET_SYSTEM_PACKAGES="libxft-dev"

	for i in $SWEET_SYSTEM_PACKAGES; do
		dpkg -s "$i" >/dev/null 2>&1
		if [ "x$?" != "x0" ]; then
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			echo ""
			echo "Ubuntu system detected and packages missing, please use"
			echo ""
			echo "    sudo apt-get install $SWEET_SYSTEM_PACKAGES"
			echo ""
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			return
		fi
	done
fi


if [ "`basename $SHELL`" != "bash" ]; then
	echo "ERROR|"
	echo "ERROR| These scripts are only compatible to the bash shell"
	echo "ERROR|"
	return
fi

export PS1="[SWEET] $PS1"

BACKDIR="$PWD"

test "x${PWD##*/}" = "xlocal_software" && cd ../

SCRIPTDIR="`pwd`/local_software"

export SWEET_ROOT="`pwd`"


if [ ! -d "$SCRIPTDIR" ]; then
	echo
	echo "ERROR| Execute this script only from the SWEET root directory"
	echo "     |   $ source local_software/env_vars.sh"
	echo
	return
fi

export PATH="$SCRIPTDIR/local/bin:$PATH"
export PKG_CONFIG_PATH="$SCRIPTDIR/local/lib/pkgconfig:$PKG_CONFIG_PATH"

export LD_LIBRARY_PATH="$SCRIPTDIR/local/lib:$LD_LIBRARY_PATH"
# TODO
# TODO check for existing directory
# TODO
export LD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$LD_LIBRARY_PATH"

export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib:$LD_LIBRARY_PATH"
export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$LD_LIBRARY_PATH"

export PYTHONPATH="$PYTHONPATH:$SCRIPTDIR/local/lib/python3.6/site-packages/"
#export PYTHONHOME="$SCRIPTDIR/local/lib/:$PYTHONHOME"

echo "SWEET environment variables loaded"

cd "$BACKDIR"

