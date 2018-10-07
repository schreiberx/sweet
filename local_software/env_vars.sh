#
# This is a template to set up environment variables correctly
#
# It assumes that you install all your libraries in subdirectories in $SWEETROOT/local_software/local
#


#######################################################################
# Backup current directory to go back to it at the end of this script
#######################################################################

BACKDIR="$PWD"


#######################################################################
# Setup echo functions ################################################
#######################################################################

#
# START: Some convenient functions
#
# echo_info [message]:
#	Output in regular colors
#
# echo_success [message]:
#	Output success message in green color
#
# echo_warning [message]:
#	Output warning message in yellow color
#
# echo_error [message]:
#	Output success message in red color
#
# echo_*_hline [message]:
#	Output a horizontal separator line
#

export SWEET_ECHO_PREFIX="echo -n \"SWEET: \""

# Pretty output
echo_info()( eval ${SWEET_ECHO_PREFIX}; echo "${@}"; )
echo_success()( eval ${SWEET_ECHO_PREFIX}; echo -en "\033[0;32m"; echo "${@}"; echo -en "\033[0m"; )
echo_warning()( eval ${SWEET_ECHO_PREFIX}; echo -en "\033[0;33m"; echo "${@}"; echo -en "\033[0m"; )
echo_error()( eval ${SWEET_ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; )

# hlines
echo_info_hline()( echo_info "*************************************************************************"; )
echo_success_hline()( echo_success "*************************************************************************"; )
echo_warning_hline()( echo_warning "*************************************************************************"; )
echo_error_hline()( echo_error "*************************************************************************"; )

# output error and exit
echo_error_exit(){ echo_error_hline; eval ${SWEET_ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; echo_error_hline; return 2>/dev/null; exit 1; }

export -f echo_info echo_success echo_warning echo_error
export -f echo_info_hline echo_success_hline echo_warning_hline echo_error_hline
export -f echo_error_exit


#######################################################################
# Check if SWEET environment are already setup ########################
#######################################################################

if [ "#$SWEET_ROOT" != "#" ]; then

	echo_warning "Environment variables already loaded (skipping)"
	return 2>/dev/null
	exit
fi



#######################################################################
# Check if this script is included (OK) or directly executed (FAILURE)
#######################################################################

#
# Include detection
#
#
SOURCED="true"


if [ "#$(basename -- $0)" = "#env_vars.sh" ]; then
	SOURCED="false"
fi

if [ "$SOURCED" != "true" ]; then
	echo_error_hline
	echo_error "THIS SCRIPT MAY NOT BE EXECUTED, BUT INCLUDED IN THE ENVIRONMENT VARIABLES!"
	echo_error_hline
	echo_error "Use e.g. "
	echo_error ""
	echo_error "   $ source ./env_vars.sh"
	echo_error ""
	echo_error "to setup the environment variables correctly"
	echo_error_hline
	return 2>/dev/null
	exit 1
fi

if [ "`basename -- "$SHELL"`" != "bash" ]; then
	echo_error_hline
	echo_error "These scripts are only compatible to the bash shell"
	echo_error_hline
	return
fi



#######################################################################
# Setup important directory environment variables #####################
#######################################################################

# Location of this script
SCRIPTDIR="$(dirname "${BASH_SOURCE[0]}")"
cd "$SCRIPTDIR"
# Get full path
SCRIPTDIR="$PWD"

# Get SWEET root directory
cd "../"
export TMP_SWEET_ROOT="$PWD"



#######################################################################
# Setup platform specific parts
#######################################################################

source $TMP_SWEET_ROOT/local_software/load_platform.sh $@ || exit 1

export SWEET_ROOT="$TMP_SWEET_ROOT"


#######################################################################
# Setup SWEET specific parts
#######################################################################
# This must be done after the platform specific parts to ensure
# that environment variables are correctly overwritten/extended with
# the SWEET variables
#######################################################################

# Back to local software
cd "$SCRIPTDIR"

#echo_info " + Setting up platform independent environment variables..."

export PATH="$SCRIPTDIR/local/bin:$PATH"
export PKG_CONFIG_PATH="$SCRIPTDIR/local/lib/pkgconfig:$PKG_CONFIG_PATH"

#
# setup SWEET_LD_LIBRARY_PATH which stored additional library paths
# This is used
#	1) to extend LD_LIBRARY_PATH in this script
#	2) to fix the mess in job schedulers (Cheyenne, e.g.) which
#	   otherwise prioritize their own system libraries such as libfftw.
#
export SWEET_LD_LIBRARY_PATH="$SCRIPTDIR/local/lib"
if [ -d "$SCRIPTDIR/local/lib64" ]; then
	export SWEET_LD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$SWEET_LD_LIBRARY_PATH"
fi
export LD_LIBRARY_PATH="$SWEET_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"


export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib:$LD_LIBRARY_PATH"
if [ -d "$SCRIPTDIR/local/lib64" ]; then
	export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$LD_LIBRARY_PATH"
fi

export PYTHONPATH="$PYTHONPATH:$SCRIPTDIR/local/lib/python3.6/site-packages/"

# Add SWEET python_mods to pythonpath
export PYTHONPATH="$PYTHONPATH:$SWEET_ROOT/python_mods"


# Back to local software
cd "$SCRIPTDIR"

PS_RELPATH="realpath  --relative-to=$SWEET_ROOT ./"

# Make prompt nice looking
#export PS1='\[\033[01;32m\][SWEET \u@$SWEET_PLATFORM_ID]\[\033[00m\] $($PS_RELPATH)\$ '
export PS1='\[\033[01;32m\][SWEET @ $SWEET_PLATFORM_ID]\[\033[00m\] $($PS_RELPATH)\$ '



#######################################################################
#######################################################################
#######################################################################

cd "$BACKDIR"

echo_success_hline
echo_success " SWEET environment setup successfully"
echo_success_hline

