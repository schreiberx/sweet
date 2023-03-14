#
# Default configuration file
#


##########
## GUIX ##
##########

## activate guix and load profile
source /applis/site/guix-start.sh

## remove guix folder from $HOME and create guix folder in local_software/local
rm $HOME/.guix-profile
if [ ! -d "$MULE_SOFTWARE_ROOT/local_software/local" ]; then
	mkdir "$MULE_SOFTWARE_ROOT/local_software/local"
fi

export GUIX_PROFILE="$MULE_SOFTWARE_ROOT/local_software/local/guix"

if [ -d "/var/guix/profiles/per-user/$USER" ]; then
  rm $GUIX_PROFILE
  ln -s /var/guix/profiles/per-user/$USER/sweet_llvm $GUIX_PROFILE
fi

guix install -p $GUIX_USER_PROFILE_DIR/sweet_llvm llvm@15 clang@15 libomp@15 mpich@3.3.2 cmake@3.25.1 gfortran-toolchain@10.3.0 lapack@3.9.0 glibc@2.32

export PATH=$GUIX_PROFILE/bin:$PATH

##############
## END GUIX ##
##############

export MULE_JOB_SCHEDULER_NUM_JOB_LIMITATION=30
export MULE_JOB_SCHEDULER_SLEEP_SECS=30

#
# Compiler environment
#
export F90=gfortran
export CC=clang
export CXX=clang++
export FC=$F90
export LD=ld
export MULE_LINK=$MULE_CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpicxx
export MULE_MPIF90=mpif90


export MULE_MPILINK=mpicxx
# If we link with mpif90, we have to add stdc++ for C++
#export MULE_MPILIBS=stdc++
export MULE_MPILIBS=gfortran


export MULE_CC_COMPILER=llvm
export MULE_CXX_COMPILER=llvm
export MULE_F90_COMPILER=gcc

