#!/usr/bin/env python

# Extensions script for diffpy.pdffit2

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

import glob
import os
import re
import sys
import warnings

from setuptools import Extension, setup

# Use this version when git data are not available, like in git zip archive.
# Update when tagging a new release.
FALLBACK_VERSION = "1.4.3"

MYDIR = os.path.dirname(os.path.abspath(__file__))

# Helper functions -----------------------------------------------------------


def get_compiler_type():
    """find compiler used for building extensions."""
    cc_arg = [a for a in sys.argv if a.startswith("--compiler=")]
    if cc_arg:
        compiler_type = cc_arg[-1].split("=", 1)[1]
    else:
        from distutils.ccompiler import new_compiler

        compiler_type = new_compiler().compiler_type
    return compiler_type


def get_gsl_config():
    """Return dictionary with paths to GSL library."""
    gslcfgpaths = [os.path.join(p, "gsl-config") for p in ([MYDIR] + os.environ["PATH"].split(os.pathsep))]
    gslcfgpaths = [p for p in gslcfgpaths if os.path.isfile(p)]
    rv = {"include_dirs": [], "library_dirs": []}
    if not gslcfgpaths:
        wmsg = "Cannot find gsl-config in {!r} nor in system PATH."
        warnings.warn(wmsg.format(MYDIR))
        return rv
    gslcfg = gslcfgpaths[0]
    with open(gslcfg) as fp:
        txt = fp.read()
    mprefix = re.search("(?m)^prefix=(.+)", txt)
    minclude = re.search(r"(?m)^[^#]*\s-I(\S+)", txt)
    mlibpath = re.search(r"(?m)^[^#]*\s-L(\S+)", txt)
    if not mprefix:
        emsg = "Cannot find 'prefix=' line in {}."
        raise RuntimeError(emsg.format(gslcfg))
    p = mprefix.group(1)
    inc = minclude.group(1) if minclude else (p + "/include")
    lib = mlibpath.group(1) if mlibpath else (p + "/lib")
    rv["include_dirs"] += [inc]
    rv["library_dirs"] += [lib]
    return rv


def get_gsl_config_win():
    """Return dictionary with paths to GSL library on Windows"""

    gsl_path = os.environ.get("GSL_PATH")
    if gsl_path:
        inc = os.path.join(gsl_path, "include")
        lib = os.path.join(gsl_path, "lib")
    else:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        inc = os.path.join(conda_prefix, "Library", "include")
        lib = os.path.join(conda_prefix, "Library", "lib")

    rv = {"include_dirs": [], "library_dirs": []}
    rv["include_dirs"] += [inc]
    rv["library_dirs"] += [lib]
    return rv


# ----------------------------------------------------------------------------

# compile and link options
define_macros = []
os_name = os.name
if os_name == "nt":
    gcfg = get_gsl_config_win()
else:
    gcfg = get_gsl_config()
include_dirs = [MYDIR] + gcfg["include_dirs"]
library_dirs = []
if sys.platform == "darwin":
    libraries = []
else:
    libraries = ["gsl"]
extra_objects = []
extra_compile_args = []
extra_link_args = []

compiler_type = get_compiler_type()
if compiler_type in ("unix", "cygwin", "mingw32"):
    extra_compile_args = ["-std=c++11", "-Wall", "-Wno-write-strings", "-O3", "-funroll-loops", "-ffast-math"]
    extra_objects += [
        os.path.join(p, "libgsl.a") for p in gcfg["library_dirs"] if os.path.isfile(os.path.join(p, "libgsl.a"))
    ]
elif compiler_type == "msvc":
    define_macros += [("_USE_MATH_DEFINES", None)]
    extra_compile_args = ["/EHs"]
    library_dirs += gcfg["library_dirs"]
# add optimization flags for other compilers if needed

# define extension arguments here
ext_kws = {
    "include_dirs": include_dirs,
    "libraries": libraries,
    "library_dirs": library_dirs,
    "define_macros": define_macros,
    "extra_compile_args": extra_compile_args,
    "extra_link_args": extra_link_args,
    "extra_objects": extra_objects,
}


# define extension here
def create_extensions():
    ext = Extension("diffpy.pdffit2.pdffit2", glob.glob("src/extensions/**/*.cc"), **ext_kws)
    return [ext]


setup_args = dict(
    ext_modules=[],
)

if __name__ == "__main__":
    setup_args["ext_modules"] = create_extensions()
    setup(**setup_args)

# End of file
