# This module is imported from top level diffpy setup.py.
# It has to define the following variables:
#     name, description, diffpy_deps, other_deps, setup_args
# Optional variables:
#     makefiles -- a list of Makefiles to be build before installation

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

# version
__id__ = "$Id$"

from distutils.core import Extension
from distutils import sysconfig
import sys
import os.path

thisfile = os.path.abspath(locals().get('__file__', 'setup_args.py'))
thisdir = os.path.dirname(thisfile)

def prependThisDir(files):
    return [os.path.join(thisdir, f) for f in files]

# name of this subpackage
name = "diffpy.pdffit2"
description = "PDFfit2 - real space structure refinement program."

# dependencies from diffpy
diffpy_deps = []

# third-party dependencies
other_deps = []

# create Extension

# define_macros
sys.path.insert(0, os.path.dirname(thisdir))
from version import __version__ as diffpy_version
sys.path.pop(0)
define_macros = [( 'VERSION', "%s" % diffpy_version )]

# compile and link options
extra_compile_args = []
extra_link_args = []
compiler = os.path.basename(sysconfig.get_config_var("CC") or "")
if compiler[:3] in ("gcc", "g++"):
    extra_compile_args = ['-O3', '-Wall', '-funroll-loops', '-ffast-math']
    extra_link_args = ['-lgsl', '-lgslcblas', '-lm']
if sys.platform == "win32":
    extra_link_args = ['libgsl.a']
# add optimization flags for other compilers later

# helper function for building with Makefile
def printDefines():
    for m, v in define_macros:
        print "'-D%s=%s'" % (m, v)
    return

pdffit2module = Extension('diffpy.pdffit2.pdffit2',
    prependThisDir([
        'pdffit2module/bindings.cc',
        'pdffit2module/misc.cc',
        'pdffit2module/pdffit2module.cc',
        'pdffit2module/pyexceptions.cc',
        'libpdffit2/Atom.cc',
        'libpdffit2/OutputStreams.cc',
        'libpdffit2/PeriodicTable.cc',
        'libpdffit2/PointsInSphere.cc',
        'libpdffit2/StringUtils.cc',
        'libpdffit2/fit.cc',
        'libpdffit2/gaussj.cc',
        'libpdffit2/metric.cc',
        'libpdffit2/nrutil.cc',
        'libpdffit2/output.cc',
        'libpdffit2/parser.cc',
        'libpdffit2/pdf.cc',
        'libpdffit2/pdffit.cc',
        'libpdffit2/pdflsmin.cc',
        'libpdffit2/scatlen.cc',
        'libpdffit2/stru.cc',
        ]),
    include_dirs = prependThisDir(['libpdffit2', 'pdffit2module', '.']),
    extra_compile_args = extra_compile_args,
    extra_link_args = extra_link_args,
    define_macros = define_macros,
)

# finally define setup_args
setup_args = {
    "name" : name,
    "description" : description,
    "packages" : [ "diffpy.pdffit2" ],
    "package_dir" : {
        "diffpy.pdffit2" : os.path.join(thisdir, "pdffit2")
        },
    "ext_modules" : [ pdffit2module ],
    "scripts" : prependThisDir([
        "applications/pdffit2",
        ]),
}

# End of file 
