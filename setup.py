#!/usr/bin/env python

# Installation script for diffpy.pdffit2

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

import sys
from setuptools import setup, find_packages
from setuptools import Extension
import fix_setuptools_chmod

# helper function
def get_compiler_type():
    """find compiler used for building extensions.
    """
    cc_arg = [a for a in sys.argv if a.startswith('--compiler=')]
    if cc_arg:
        compiler_type = cc_arg[-1].split('=', 1)[1]
    else:
        from distutils.ccompiler import new_compiler
        compiler_type = new_compiler().compiler_type
    return compiler_type

# compile and link options
extra_compile_args = []
extra_link_args = []
libraries = ['gsl']

compiler_type = get_compiler_type()
if compiler_type in ("unix", "cygwin", "mingw32"):
    extra_compile_args = [ '-Wall', '-Wno-write-strings',
            '-O3', '-funroll-loops', '-ffast-math' ]
    libraries += ['gslcblas', 'm']
elif compiler_type == "msvc":
    extra_compile_args = [ '/EHs' ]
# add optimization flags for other compilers later

# define extension here
pdffit2module = Extension('diffpy.pdffit2.pdffit2', [
            'pdffit2module/bindings.cc',
            'pdffit2module/misc.cc',
            'pdffit2module/pdffit2module.cc',
            'pdffit2module/pyexceptions.cc',
            'libpdffit2/Atom.cc',
            'libpdffit2/LocalPeriodicTable.cc',
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
            ],
        include_dirs = ['libpdffit2', 'pdffit2module', '.'],
        extra_compile_args = extra_compile_args,
        extra_link_args = extra_link_args,
        libraries = libraries,
)

# figure which scripts need to be installed
scripts = ['applications/pdffit2']
if sys.platform.startswith('win32'):
    scripts.append('applications/pdffit2.bat')

# define distribution
setup(
        name = 'diffpy.pdffit2',
        version = '1.0',
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        scripts = scripts,
        ext_modules = [pdffit2module],
        install_requires = [
            'diffpy.Structure>=1.1-r6772',
        ],
        dependency_links = [
            'http://www.diffpy.org/packages/',
        ],

        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        maintainer = 'Pavol Juhas',
        maintainer_email = 'pj2192@columbia.edu',
        url = 'http://www.diffpy.org/',
        download_url = 'http://www.diffpy.org/packages/',
        description = 'PDFfit2 - real space structure refinement program.',
        license = 'BSD',
        keywords = 'PDF structure refinement',
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: C++',
            'Programming Language :: Python :: 2.5',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file
