#!/usr/bin/env python

# This script is usually called from top level diffpy installation.

"""PDFfit2 - real space structure refinement engine

Installs to:    diffpy.pdffit2
"""

# version
__id__ = "$Id$"

from distutils.core import setup
from distutils.core import Extension
from distutils import sysconfig
import sys
import os

thisfile = os.path.abspath(locals().get('__file__', 'setup.py'))
setup_dir = os.path.dirname(thisfile)

sys.path.insert(0, setup_dir)
import pdffit2.version
sys.path.pop(0)
package_version = pdffit2.version.__version__
package_date = pdffit2.version.__date__

define_macros = [( 'VERSION', '"%s"' % package_version )]

# figure out optimization options
extra_compile_args = []
compiler = os.path.basename(sysconfig.get_config_var("CC"))
if compiler[:3] in ("gcc", "g++"):
    extra_compile_args = ['-O3', '-Wall', '-funroll-loops', '-ffast-math']
# add optimization flags for other compilers later

def printDefines():
    for m, v in define_macros:
        print "'-D%s=%s'" % (m, v)
    return

def prependSetupDir(files):
    return [os.path.join(setup_dir, f) for f in files]

pdffit2module = Extension('pdffit2module',
    prependSetupDir([
        'pdffit2module/bindings.cc',
        'pdffit2module/exceptions.cc',
        'pdffit2module/misc.cc',
        'pdffit2module/pdffit2module.cc',
        'libpdffit2/Atom.cc',
        'libpdffit2/PeriodicTable.cc',
        'libpdffit2/PointsInSphere.cc',
        'libpdffit2/StringUtils.cc',
        'libpdffit2/fit.cc',
        'libpdffit2/gaussj.cc',
        'libpdffit2/math.cc',
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
    include_dirs = prependSetupDir(['libpdffit2', 'pdffit2module', '.']),
    extra_compile_args = extra_compile_args,
    extra_link_args = ['-lgsl', '-lgslcblas', '-lm'],
    define_macros = define_macros,
)

# define distribution
setup_args = {
    "name" : "pdffit2",
    "description" : "PDFfit2 - real space structure refinement program.",
    "version" : package_version,
    "packages" : [ "diffpy.pdffit2" ],
    "package_dir" : {"diffpy.pdffit2" : os.path.join(setup_dir, "pdffit2")},
    "ext_package" : 'diffpy.pdffit2',
    "ext_modules" : [ pdffit2module ],
    "scripts" : prependSetupDir(["applications/pdffit2"]),
}

diffpy__init__code = """
import sys
import os.path
mydir = os.path.dirname(__file__)
if mydir not in sys.path:   sys.path.insert(0, mydir)
""".lstrip()

def check_diffpy__init__(distribution):
    """check if diffpy has __init__.py and create one if not
    """
    from distutils import log
    install_lib = None
    if 'install' in distribution.commands:
        opts = distribution.get_option_dict('install')
        install_lib = opts.get('install_lib', 2*[None])[1]
    if 'install_lib' in distribution.commands and not install_lib:
        opts = distribution.get_option_dict('install_lib')
        install_lib = opts.get('install_dir', 2*[None])[1]
    if not install_lib:             return
    initfile = os.path.join(install_lib, 'diffpy', '__init__.py')
    if os.path.isfile(initfile):    return
    # we need to create and compile the file
    log.info("creating " + initfile)
    out = open(initfile, 'w')
    out.write(diffpy__init__code)
    out.close()
    import compiler
    log.info("byte-compiling %s to %s" % \
            (initfile, os.path.basename(initfile)) )
    compiler.compileFile(initfile)
    return

if __name__ == "__main__":
    distribution = setup(**setup_args)
    check_diffpy__init__(distribution)

# End of file 
