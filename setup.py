#!/usr/bin/env python

# Extensions script for diffpy.pdffit2

import glob
import sys

from setuptools import Extension, setup

# Define extension arguments here
ext_kws = {
    "libraries": [],
    "extra_compile_args": [],
    "extra_link_args": [],
    "include_dirs": [],
}


# Figure out the tagged name of boost_python library.
def get_boost_libraries():
    """Check for installed boost_python shared library.

    Returns list of required boost_python shared libraries that are installed
    on the system. If required libraries are not found, an Exception will be
    thrown.
    """
    baselib = "boost_python"
    major, minor = (str(x) for x in sys.version_info[:2])
    pytags = [major + minor, major, ""]
    mttags = ["", "-mt"]
    boostlibtags = [(pt + mt) for mt in mttags for pt in pytags] + [""]
    from ctypes.util import find_library

    for tag in boostlibtags:
        lib = baselib + tag
        found = find_library(lib)
        if found:
            break

    # Show warning when library was not detected.
    if not found:
        import platform
        import warnings

        ldevname = "LIBRARY_PATH"
        if platform.system() == "Darwin":
            ldevname = "DYLD_FALLBACK_LIBRARY_PATH"
        wmsg = ("Cannot detect name suffix for the %r library. " "Consider setting %s.") % (baselib, ldevname)
        warnings.warn(wmsg)

    libs = [lib]
    return libs


def create_extensions():
    "Initialize Extension objects for the setup function."
    blibs = [n for n in get_boost_libraries() if n not in ext_kws["libraries"]]
    ext_kws["libraries"] += blibs
    ext = Extension("diffpy.pdffit2.pdffit2_ext", glob.glob("src/extensions/*.cpp"), **ext_kws)
    return [ext]


# Extensions not included in pyproject.toml
setup_args = dict(
    ext_modules=[],
)


if __name__ == "__main__":
    setup_args["ext_modules"] = create_extensions()
    setup(**setup_args)
