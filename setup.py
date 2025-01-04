#!/usr/bin/env python

# Extensions script for diffpy.pdffit2

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

import glob
import os
import re
import shutil
import sys
import warnings
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# Use this version when git data are not available, like in git zip archive.
# Update when tagging a new release.
FALLBACK_VERSION = "1.4.3"

MYDIR = str(Path(__file__).parent.resolve())

# Helper functions -----------------------------------------------------------


def get_compiler_type():
    """Find compiler used for building extensions."""
    cc_arg = [a for a in sys.argv if a.startswith("--compiler=")]
    if cc_arg:
        return cc_arg[-1].split("=", 1)[1]
    from distutils.ccompiler import new_compiler

    return new_compiler().compiler_type


def get_gsl_config():
    """Return dictionary with paths to GSL library."""
    gslcfgpaths = [Path(p) / "gsl-config" for p in ([MYDIR] + os.environ["PATH"].split(os.pathsep))]
    gslcfgpaths = [p for p in gslcfgpaths if p.is_file()]
    rv = {"include_dirs": [], "library_dirs": []}
    if not gslcfgpaths:
        warnings.warn(f"Cannot find gsl-config in {MYDIR} nor in system PATH.")
        return rv
    gslcfg = gslcfgpaths[0]
    txt = gslcfg.read_text()
    mprefix = re.search(r"(?m)^prefix=(.+)", txt)
    minclude = re.search(r"(?m)^[^#]*\s-I(\S+)", txt)
    mlibpath = re.search(r"(?m)^[^#]*\s-L(\S+)", txt)
    if not mprefix:
        raise RuntimeError(f"Cannot find 'prefix=' line in {gslcfg}.")
    p = Path(mprefix.group(1))
    rv["include_dirs"].append(str(minclude.group(1) if minclude else p / "include"))
    rv["library_dirs"].append(str(mlibpath.group(1) if mlibpath else p / "lib"))
    return rv


def get_gsl_config_win():
    """Return dictionary with paths to GSL library on Windows."""
    gsl_path = os.environ.get("GSL_PATH", "")
    if gsl_path:
        inc = Path(gsl_path) / "include"
        lib = Path(gsl_path) / "lib"
    else:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            inc = Path(conda_prefix) / "Library" / "include"
            lib = Path(conda_prefix) / "Library" / "lib"
        else:
            raise EnvironmentError(
                "Neither GSL_PATH nor CONDA_PREFIX environment variables are set. "
                "Please ensure GSL is installed and GSL_PATH is correctly set."
            )
    return {"include_dirs": [str(inc)], "library_dirs": [str(lib)]}


class CustomBuildExt(build_ext):
    def run(self):
        super().run()
        gsl_path = (
            Path(os.environ.get("GSL_PATH"))
            if os.environ.get("GSL_PATH")
            else Path(os.environ.get("CONDA_PREFIX", "")) / "Library"
        )
        bin_path = gsl_path / "bin"
        dest_path = Path(self.build_lib) / "diffpy" / "pdffit2"
        dest_path.mkdir(parents=True, exist_ok=True)

        for dll_file in bin_path.glob("gsl*.dll"):
            shutil.copy(str(dll_file), str(dest_path))


# ----------------------------------------------------------------------------

# Compile and link options
os_name = os.name
if os_name == "nt":
    gcfg = get_gsl_config_win()
else:
    gcfg = get_gsl_config()

if sys.platform == "darwin":
    libraries = []
else:
    libraries = ["gsl"]

include_dirs = [MYDIR] + gcfg["include_dirs"]
library_dirs = []
define_macros = []
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

# Define extension arguments
ext_kws = {
    "include_dirs": include_dirs,
    "libraries": libraries,
    "library_dirs": library_dirs,
    "define_macros": define_macros,
    "extra_compile_args": extra_compile_args,
    "extra_link_args": extra_link_args,
    "extra_objects": extra_objects,
}


# Define extensions
def create_extensions():
    ext = Extension("diffpy.pdffit2.pdffit2", glob.glob("src/extensions/**/*.cc"), **ext_kws)
    return [ext]


setup_args = dict(
    ext_modules=[],
    cmdclass={"build_ext": CustomBuildExt},
)

if __name__ == "__main__":
    setup_args["ext_modules"] = create_extensions()
    setup(**setup_args)

# End of file
