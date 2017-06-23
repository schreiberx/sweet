#! /bin/bash

EXEC="gdb -ex run -ex bt --args $@"

echo "$EXEC"
$EXEC
