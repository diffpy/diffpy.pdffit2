#!/bin/bash

# Mac OS X Python 2.6 needs explicit link with the standard C++ library.
if [[ "$(uname)" == Darwin && "$PY_VER" == 2.6 ]]; then
    LDFLAGS="${LDFLAGS} -lstdc++"
fi

# Produce stripped binaries on Linux.
if [[ "$(uname)" == Linux ]]; then
    # conda-build does not set LDFLAGS on Linux.
    export LDFLAGS="-s"
fi

$PYTHON setup.py install

# Add more build steps here, if they are necessary.

$PYTHON setup.py --version > __conda_version__.txt

# See http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
