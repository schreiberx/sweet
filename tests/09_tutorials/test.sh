#! /bin/bash

set -e

cd "$MULE_SOFTWARE_ROOT"

cd tutorials

./run_all_tutorials.sh

./run_cleanup.sh

