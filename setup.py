#!/usr/bin/env python

# Extensions script for diffpy.pdffit2

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

import os
import re
import sys
import warnings

from setuptools import Extension, find_packages, setup

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
    """Return dictionary with paths to GSL library, windwows version.
    This version is installed with conda.
    """
    conda_prefix = os.environ["CONDA_PREFIX"]
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
libraries = []
extra_objects = []
extra_compile_args = []
extra_link_args = []

compiler_type = get_compiler_type()
if compiler_type in ("unix", "cygwin", "mingw32"):
    extra_compile_args = ["-std=c++11", "-Wall", "-Wno-write-strings", "-O3", "-funroll-loops", "-ffast-math"]
    extra_objects += ((p + "/libgsl.a") for p in gcfg["library_dirs"])
elif compiler_type == "msvc":
    define_macros += [("_USE_MATH_DEFINES", None)]
    extra_compile_args = ["/EHs"]
    libraries += ["gsl"]
    library_dirs += gcfg["library_dirs"]
# add optimization flags for other compilers if needed


# define extension here
pdffit2module = Extension(
    "diffpy.pdffit2.pdffit2",
    [
        "src/extensions/pdffit2module/bindings.cc",
        "src/extensions/pdffit2module/misc.cc",
        "src/extensions/pdffit2module/pdffit2module.cc",
        "src/extensions/pdffit2module/pyexceptions.cc",
        "src/extensions/libpdffit2/Atom.cc",
        "src/extensions/libpdffit2/LocalPeriodicTable.cc",
        "src/extensions/libpdffit2/OutputStreams.cc",
        "src/extensions/libpdffit2/PeriodicTable.cc",
        "src/extensions/libpdffit2/PointsInSphere.cc",
        "src/extensions/libpdffit2/StringUtils.cc",
        "src/extensions/libpdffit2/fit.cc",
        "src/extensions/libpdffit2/gaussj.cc",
        "src/extensions/libpdffit2/metric.cc",
        "src/extensions/libpdffit2/nrutil.cc",
        "src/extensions/libpdffit2/output.cc",
        "src/extensions/libpdffit2/parser.cc",
        "src/extensions/libpdffit2/pdf.cc",
        "src/extensions/libpdffit2/pdffit.cc",
        "src/extensions/libpdffit2/pdflsmin.cc",
        "src/extensions/libpdffit2/scatlen.cc",
        "src/extensions/libpdffit2/stru.cc",
    ],
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs,
    define_macros=define_macros,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    extra_objects=extra_objects,
)

setup_args = dict(
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=[pdffit2module],
    # scripts=[]    # place examples here
)

if __name__ == "__main__":
    setup(**setup_args)

# End of file
