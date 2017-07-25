#!/bin/bash

# Produce stripped binaries on Linux.
if [[ "$(uname)" == Linux ]]; then
    # conda-build does not set LDFLAGS on Linux.
    export LDFLAGS="-s"
fi

$PYTHON setup.py install

# Add more build steps here, if they are necessary.

# See http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
