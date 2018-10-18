#! /bin/bash

cd "$SWEET_ROOT"

SCONS="scons --program=burgers --parareal=serial --plane-spectral-space=enable --mode=debug"
echo $SCONS
$SCONS || exit

