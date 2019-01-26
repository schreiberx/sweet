#! /usr/bin/env bash

BENCHMARK_DIRNAME=$1

if [[ -z "$BENCHMARK_DIRNAME" ]]; then

	if [[ -e "benchmarks_create.py" ]]; then
		# We are currently in the benchmark directory itself
		CURDIR="`pwd`"
		BENCHMARK_DIRNAME=`basename "$CURDIR"`

		echo ""
		echo "Using current directory '$BENCHMARK_DIRNAME' for benchmark directory"
		echo ""

		echo "Changing current directory to '../'"
		echo ""
		cd ..

	else
		echo ""
		echo "Create a tarball of one benchmkark directory (with job directories as subfolders)"
		echo "All large output files (*.csv) will be excluded"
		echo ""
		echo "Usage:"
		echo "	$0 [benchmark_folder_name]"
		echo "	where 'benchmark_folder_name' is the folder in the current directory to be tarballed."
		echo ""
		echo "or"
		echo "	$0"
		echo "  which is to be executed within the benchmark folder to be tarballed."
		echo "  a tarball is then created in the upper directory."
		echo ""

	fi
fi


# Full path to benchmark directory
#BENCHMARK_PATH="$(realpath $BENCHMARK_DIRNAME)"

# Ending of TARBALL
TARENDING="tar.xz"

# Tar options
TAROPTIONS="cJvf"

# Full path to tarball
#TARBALL_PATH="$BENCHMARK_PATH.$TARENDING"

# Filename of tarball
TARFILENAME="$BENCHMARK_DIRNAME.$TARENDING"


if [[ -e "$TARBALL" ]]; then
	echo "Error: File $TARBALL already exists"
	exit 1
fi

echo "Creating tarball $TARBALL"
EXEC="tar $TAROPTIONS $TARFILENAME $BENCHMARK_DIRNAME --exclude=*csv"
echo " \$ $EXEC"

$EXEC

