
#
# Automatically detect the most recent clang++ version
#
# Start at version 30 and search downwards
#


type "clang++" 2> /dev/null 1>&2
if [[ "$?" == "0" ]]; then
	CLANG_POSTFIX=""
else

	CLANG_VERSIONS=$(seq 30 -1 0)
	for i in $CLANG_VERSIONS; do
		type "clang++-$i" 2> /dev/null 1>&2

		if [[ "$?" == "0" ]]; then
			CLANG_VERSION=$i
			CLANG_POSTFIX="-$CLANG_VERSION"
			break
		fi
	done

	if [[ "$i" == "0" ]]; then
		CLANG_POSTFIX=""
		#echo "No clang++ compiler found!"
		#echo ""
		#echo "Not setting up any variables."
		#echo ""
		#exit 0
	fi

fi

test -z "$CC" && export CC=clang$CLANG_POSTFIX
test -z "$CXX" && export CXX=clang++$CLANG_POSTFIX

# We still use gfortran
test -z "$F90" && export F90=gfortran
test -z "$FC" && export FC=$F90

# Use standard LD
test -z "$LD" && export LD=ld


# Setup further MULE variables
test -z "$MULE_LINK" && export MULE_LINK=$MULE_CXX

test -z "$MULE_MPICC" && export MULE_MPICC=mpicc
test -z "$MULE_MPICXX" && export MULE_MPICXX=mpic++
test -z "$MULE_MPI90" && export MULE_MPIF90=mpif90

test -z "$MULE_MPILINK" && export MULE_MPILINK=mpif90

# If we link with mpif90, we have to add stdc++ for C++
#test -z "MULE_MPILIBS" && export MULE_MPILIBS=stdc++


test -z "$MULE_CC_COMPILER" && export MULE_CC_COMPILER="llvm"
test -z "$MULE_CXX_COMPILER" && export MULE_CXX_COMPILER="llvm"
test -z "$MULE_F90_COMPILER" && export MULE_F90_COMPILER="gcc"

