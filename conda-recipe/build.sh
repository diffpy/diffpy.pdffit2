#!/bin/bash

# Produce stripped binaries on Linux.
if [[ "$(uname)" == Linux ]]; then
    # conda-build does not set LDFLAGS on Linux.
    export LDFLAGS="-s"
fi

# allow use of default SDK on travis
if [[ "$(id -un)" == travis ]]; then
    unset CONDA_BUILD_SYSROOT
# print informative message if SDK is not installed
elif [[ "$(uname)" == Darwin ]]; then
    if [[ ! -e $CONDA_BUILD_SYSROOT ]]; then
	echo "Please install macOS SDK to $CONDA_BUILD_SYSROOT"
	exit 2
    fi
fi

$PYTHON -m easy_install --no-deps .

# Add more build steps here, if they are necessary.

# See http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
