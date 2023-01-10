#! /bin/bash

source activate.sh || exit 1

cd tests

./test.sh
