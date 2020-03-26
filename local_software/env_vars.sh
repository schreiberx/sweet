#
# This is a template to set up environment variables correctly
#
# It assumes that you install all your libraries in subdirectories in $MULE_SOFTWARE_ROOT/local_software/local
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

export MULE_ECHO_PREFIX="echo -n \"MULE: \""

# Pretty output
echo_info()( eval ${MULE_ECHO_PREFIX}; echo "${@}"; )
echo_success()( eval ${MULE_ECHO_PREFIX}; echo -en "\033[0;32m"; echo "${@}"; echo -en "\033[0m"; )
echo_warning()( eval ${MULE_ECHO_PREFIX}; echo -en "\033[0;33m"; echo "${@}"; echo -en "\033[0m"; )
echo_error()( eval ${MULE_ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; )

# hlines
echo_info_hline()( echo_info "*************************************************************************"; )
echo_success_hline()( echo_success "*************************************************************************"; )
echo_warning_hline()( echo_warning "*************************************************************************"; )
echo_error_hline()( echo_error "*************************************************************************"; )

# output error and exit
echo_error_exit(){ echo_error_hline; eval ${MULE_ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; echo_error_hline; exit 1; }
#echo_error_return(){ echo_error_hline; eval ${MULE_ECHO_PREFIX}; echo -en "\033[0;31m"; echo "${@}"; echo -en "\033[0m"; echo_error_hline; return; }

echo_exec()( eval ${MULE_ECHO_PREFIX}; echo "Executing '${@}'"; ${@})

export -f echo_info echo_success echo_warning echo_error
export -f echo_info_hline echo_success_hline echo_warning_hline echo_error_hline
export -f echo_error_exit
export -f echo_exec


#######################################################################
# Check if MULE_SOFTWARE_ROOT environment are already setup ########################
#######################################################################

if [ "#$MULE_SOFTWARE_ROOT" != "#" ]; then

	echo_warning "Environment variables already loaded for platform '${MULE_PLATFORM_ID}'(skipping)"
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

# Get SOFTWARE root directory
cd "../"
export MULE_SOFTWARE_ROOT="$PWD"

# Setup mule main path
export MULE_ROOT="$MULE_SOFTWARE_ROOT/mule"
export MULE_LOCAL_ROOT="$MULE_SOFTWARE_ROOT/mule_local"


#######################################################################
# Setup platform specific parts
#######################################################################

source $MULE_ROOT/bin/load_platform.sh $@ || exit 1

if [[ -z "$MULE_PLATFORM_DIR" ]]; then
	unset MULE_ROOT
	unset MULE_LOCAL_ROOT
	unset MULE_SOFTWARE_ROOT
	unset SCRIPTDIR
	return
fi


#######################################################################
# Setup SOFTWARE specific parts
#######################################################################
# This must be done after the platform specific parts to ensure
# that environment variables are correctly overwritten/extended with
# the SWEET variables
#######################################################################

# Back to local software
cd "$SCRIPTDIR"

#echo_info " + Setting up platform independent environment variables..."

export PATH="$SCRIPTDIR/local/bin:$PATH"

export PATH="$MULE_ROOT/bin:$PATH"
export PATH="$MULE_LOCAL_ROOT/bin:$PATH"
export PATH="$MULE_SOFTWARE_ROOT/bin:$PATH"

export PKG_CONFIG_PATH="$SCRIPTDIR/local/lib/pkgconfig:$PKG_CONFIG_PATH"


#
# setup MULE_LD_LIBRARY_PATH which stored additional library paths
# This is used
#	1) to extend LD_LIBRARY_PATH in this script
#	2) to fix the mess in job schedulers (Cheyenne, e.g.) which
#	   otherwise prioritize their own system libraries such as libfftw.
#
export MULE_LD_LIBRARY_PATH="$SCRIPTDIR/local/lib"
# Always include lib64 even if it doesn exist.
# Thats important to install software into this directory
export MULE_LD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$MULE_LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$MULE_LD_LIBRARY_PATH:$LD_LIBRARY_PATH"


export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib:$LD_LIBRARY_PATH"
if [ -d "$SCRIPTDIR/local/lib64" ]; then
	export DYLD_LIBRARY_PATH="$SCRIPTDIR/local/lib64:$LD_LIBRARY_PATH"
fi

PYTHONVERSION=$(python3 -c "import sys;print(str(sys.version_info.major)+\".\"+str(sys.version_info.minor),end='')")
export PYTHONPATH="$SCRIPTDIR/local/lib/python$PYTHONVERSION/site-packages/:$PYTHONPATH"


# Add MULE python path
# DO NOT INCLUDE THIS, SINCE ONLY ONE 'mule' MODULE WILL BE AVAILABLE!
# We link from the SWEET module path to the MULE module path
export PYTHONPATH="$PYTHONPATH:$MULE_ROOT/python"

# Add MULE SWEET python path
export PYTHONPATH="$PYTHONPATH:$MULE_LOCAL_ROOT/python"


# Back to local software
cd "$SCRIPTDIR"

# Test if 'realpath' exists
type realpath >/dev/null 2>&1
if [[ $? -eq 0 ]]; then
	SWEET_SHELL_PATH='$(realpath  --relative-base=$MULE_SOFTWARE_ROOT ./)'
	#export PS1='\[\033[01;32m\][SWEET \u@$MULE_PLATFORM_ID]\[\033[00m\] $($PS_RELPATH)\$ '
	export PS1="\[\033[01;32m\][SWEET @ $MULE_PLATFORM_ID]\[\033[00m\] $SWEET_SHELL_PATH\$ "
else
	SWEET_SHELL_PATH='\w'
	export PS1="\[\033[01;32m\][SWEET @ $MULE_PLATFORM_ID]\[\033[00m\] $SWEET_SHELL_PATH\$ "
fi



#######################################################################
#######################################################################
#######################################################################

cd "$BACKDIR"

echo_success_hline
echo_success " MULE SOFTWARE environment setup successfully"
echo_success_hline

