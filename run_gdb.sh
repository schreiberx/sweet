#! /bin/bash

EXEC="gdb -d ./ -ex run -ex bt --args $@"

echo "$EXEC"
$EXEC
