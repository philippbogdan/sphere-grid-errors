#!/bin/sh
# Compile the checked NCAR source with the local grid writer and 64-bit reals.
set -eu
if [ "$#" -ne 4 ]; then
    echo "Usage: sh scripts/build.sh SOURCE_DIRECTORY TRUNCATION NLAT NLON" >&2
    exit 2
fi
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
"${PYTHON:-python3}" "$SCRIPT_DIR/build_reference.py" "$1" "$2" "$3" "$4"
