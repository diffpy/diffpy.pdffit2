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
    """Return the compiler type used during the build."""
    cc_arg = [a for a in sys.argv if a.startswith("--compiler=")]
    if cc_arg:
        return cc_arg[-1].split("=", 1)[1]
    from distutils.ccompiler import new_compiler

    return new_compiler().compiler_type


def get_gsl_config():
    """
    Determine the GSL include and library directories by trying in order:
      1) CONDA_PREFIX,
      2) GSL_PATH,
      3) gsl-config (for Unix-like systems).
    Raises EnvironmentError if none are found.
    """
    rv = {"include_dirs": [], "library_dirs": []}

    # 1. Check using CONDA_PREFIX.
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        if os.name == "nt":
            inc = Path(conda_prefix) / "Library" / "include"
            lib = Path(conda_prefix) / "Library" / "lib"
        else:
            inc = Path(conda_prefix) / "include"
            lib = Path(conda_prefix) / "lib"
        if inc.is_dir() and lib.is_dir():
            rv["include_dirs"].append(str(inc))
            rv["library_dirs"].append(str(lib))
            return rv
        else:
            warnings.warn(
                f"CONDA_PREFIX is set to {conda_prefix}, " "but GSL not found at those paths. Proceeding..."
            )

    # 2. Check using GSL_PATH.
    gsl_path = os.environ.get("GSL_PATH", "")
    if gsl_path:
        inc = Path(gsl_path) / "include"
        lib = Path(gsl_path) / "lib"
        if inc.is_dir() and lib.is_dir():
            rv["include_dirs"].append(str(inc))
            rv["library_dirs"].append(str(lib))
            return rv
        else:
            raise EnvironmentError(
                f"GSL_PATH={gsl_path} is set, but {inc} or {lib} not found. " "Please verify your GSL_PATH."
            )

    # 3. Try using the gsl-config executable (only on Unix-like systems).
    if os.name != "nt":
        path_dirs = os.environ.get("PATH", "").split(os.pathsep)
        gslcfg_paths = [Path(p) / "gsl-config" for p in path_dirs if p]
        gslcfg_paths = [p for p in gslcfg_paths if p.is_file()]
        if gslcfg_paths:
            gslcfg = gslcfg_paths[0]
            txt = gslcfg.read_text()
            prefix_match = re.search(r"(?m)^prefix=(.+)", txt)
            include_match = re.search(r"(?m)^[^#]*\s-I(\S+)", txt)
            lib_match = re.search(r"(?m)^[^#]*\s-L(\S+)", txt)
            if prefix_match:
                prefix_path = Path(prefix_match.group(1))
                inc_dir = include_match.group(1) if include_match else (prefix_path / "include")
                lib_dir = lib_match.group(1) if lib_match else (prefix_path / "lib")
                rv["include_dirs"].append(str(inc_dir))
                rv["library_dirs"].append(str(lib_dir))
                return rv
            else:
                raise RuntimeError(f"Cannot parse 'prefix=' from {gslcfg}.")
        else:
            warnings.warn(
                "No gsl-config found in PATH. GSL may not be installed or not in PATH. "
                "Proceeding without GSL configuration."
            )

    # 4. Nothing found: raise error.
    raise EnvironmentError(
        "Unable to locate GSL:\n"
        "1) CONDA_PREFIX not set or no GSL there\n"
        "2) GSL_PATH not set or invalid\n"
        "3) gsl-config not available\n"
        "Please set GSL_PATH or use a conda environment with GSL."
    )


class CustomBuildExt(build_ext):
    def run(self):
        # Retrieve the GSL library directories and append them to each extension.
        gsl_cfg = get_gsl_config()
        lib_dirs = gsl_cfg.get("library_dirs", [])
        for ext in self.extensions:
            # Add gsl lib for linking.
            ext.library_dirs.extend(lib_dirs)
            # Embed RPATH flags, runtime linking without LD_LIBRARY_PATH.
            ext.extra_link_args = ext.extra_link_args or []
            for lib in lib_dirs:
                ext.extra_link_args.append(f"-Wl,-rpath,{lib}")
        super().run()
        # Avoid dll error
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


def create_extensions():
    """Create the list of Extension objects for the build."""
    # lazy evaluation prevents build sdist failure
    try:
        gcfg = get_gsl_config()
    except EnvironmentError:
        return []

    # On macOS, dynamic linking may not be needed
    if sys.platform == "darwin":
        libraries = []
    else:
        libraries = ["gsl"]

    include_dirs = [MYDIR] + gcfg["include_dirs"]
    library_dirs = gcfg["library_dirs"]
    define_macros = []
    extra_objects = []
    extra_compile_args = []
    extra_link_args = []

    compiler_type = get_compiler_type()
    if compiler_type in ("unix", "cygwin", "mingw32"):
        extra_compile_args = ["-std=c++11", "-Wall", "-Wno-write-strings", "-O3", "-funroll-loops", "-ffast-math"]
        # Check for static GSL libraries and add them if found.
        static_libs = [
            os.path.join(p, "libgsl.a")
            for p in gcfg["library_dirs"]
            if os.path.isfile(os.path.join(p, "libgsl.a"))
        ]
        if static_libs:
            extra_objects += static_libs
            # Use static linking: remove "-lgsl" to avoid dynamic linking conflicts.
            libraries = []
    elif compiler_type == "msvc":
        define_macros += [("_USE_MATH_DEFINES", None)]
        extra_compile_args = ["/EHs"]

    # Extension keyword arguments.
    ext_kws = {
        "include_dirs": include_dirs,
        "libraries": libraries,
        "library_dirs": library_dirs,
        "define_macros": define_macros,
        "extra_compile_args": extra_compile_args,
        "extra_link_args": extra_link_args,
        "extra_objects": extra_objects,
    }
    ext = Extension("diffpy.pdffit2.pdffit2", glob.glob("src/extensions/**/*.cc"), **ext_kws)
    return [ext]


setup_args = dict(
    ext_modules=[],
    cmdclass={"build_ext": CustomBuildExt},
)

if __name__ == "__main__":
    setup_args["ext_modules"] = create_extensions()
    setup(**setup_args)
