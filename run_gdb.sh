#! /bin/bash

EXEC="gdb -d ./"

# Directly run debugger
EXEC+=" -ex run"

# Get backtrace
EXEC+=" -ex backtrace"

# Deactivate quit confirmation
EXEC+=" -ex stop"


# Program arguments
EXEC+=" --args $@"

echo "$EXEC"
$EXEC
