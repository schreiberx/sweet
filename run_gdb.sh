#! /bin/bash

EXEC="gdb -d ./"

# Directly run debugger
EXEC+=" -ex run"

# Get backtrace
EXEC+=" -ex bt"

# Quit debugger if program finishes
EXEC+=" -ex q"

# Program arguments
EXEC+=" --args $@"

echo "$EXEC"
$EXEC
