#!/usr/bin/env python

# Installation script for diffpy.pdffit2

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

import os
import re
import sys
import warnings

from setuptools import setup, find_packages
from setuptools import Extension

# Use this version when git data are not available, like in git zip archive.
# Update when tagging a new release.
FALLBACK_VERSION = '1.4.1'

# determine if we run with Python 3.
PY3 = (sys.version_info[0] == 3)

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
MYDIR = os.path.dirname(os.path.abspath(__file__))
versioncfgfile = os.path.join(MYDIR, 'diffpy/pdffit2/version.cfg')
gitarchivecfgfile = os.path.join(MYDIR, '.gitarchive.cfg')


def gitinfo():
    from subprocess import Popen, PIPE, check_output
    kw = dict(stdout=PIPE, cwd=MYDIR, universal_newlines=True)
    proc = Popen(['git', 'describe', '--tags', '--match=[v,V,[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %ct %ci'], **kw)
    glog = proc.stdout.read()
    rv = {}
    rv['commit'], rv['timestamp'], rv['date'] = glog.strip().split(None, 2)
    version = '.post'.join(desc.strip().split('-')[:2]).lstrip('vV')
    rv['version'] = version
    return rv


def getversioncfg():
    if PY3:
        from configparser import RawConfigParser
    else:
        from ConfigParser import RawConfigParser
    vd0 = dict(version=FALLBACK_VERSION, commit='', date='', timestamp=0)
    # first fetch data from gitarchivecfgfile, ignore if it is unexpanded
    g = vd0.copy()
    cp0 = RawConfigParser(vd0)
    cp0.read(gitarchivecfgfile)
    if len(cp0.get('DEFAULT', 'commit')) > 20:
        g = cp0.defaults()
        mx = re.search(r'\btag: [vV]?(\d[^,]*)', g.pop('refnames'))
        if mx:
            g['version'] = mx.group(1)
    # then try to obtain version data from git.
    gitdir = os.path.join(MYDIR, '.git')
    if os.path.exists(gitdir) or 'GIT_DIR' in os.environ:
        try:
            g = gitinfo()
        except OSError:
            pass
    # finally, check and update the active version file

    cp = RawConfigParser()
    cp.read(versioncfgfile)
    d = cp.defaults()
    rewrite = not d or (g['commit'] and (
        g['version'] != d.get('version') or g['commit'] != d.get('commit')))
    if rewrite:
        cp.set('DEFAULT', 'version', g['version'])
        cp.set('DEFAULT', 'commit', g['commit'])
        cp.set('DEFAULT', 'date', g['date'])
        cp.set('DEFAULT', 'timestamp', g['timestamp'])
        with open(versioncfgfile, 'w') as fp:
            cp.write(fp)
    return cp

versiondata = getversioncfg()

# Helper functions -----------------------------------------------------------

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


def get_gsl_config():
    '''Return dictionary with paths to GSL library.
    '''
    gslcfgpaths = [os.path.join(p, 'gsl-config')
                   for p in ([MYDIR] + os.environ['PATH'].split(os.pathsep))]
    gslcfgpaths = [p for p in gslcfgpaths if os.path.isfile(p)]
    rv = {'include_dirs': [], 'library_dirs': []}
    if not gslcfgpaths:
        wmsg = "Cannot find gsl-config in {!r} nor in system PATH."
        warnings.warn(wmsg.format(MYDIR))
        return rv
    gslcfg = gslcfgpaths[0]
    with open(gslcfg) as fp:
        txt = fp.read()
    mprefix = re.search('(?m)^prefix=(.+)', txt)
    minclude = re.search(r'(?m)^[^#]*\s-I(\S+)', txt)
    mlibpath = re.search(r'(?m)^[^#]*\s-L(\S+)', txt)
    if not mprefix:
        emsg = "Cannot find 'prefix=' line in {}."
        raise RuntimeError(emsg.format(gslcfg))
    p = mprefix.group(1)
    inc = minclude.group(1) if minclude else (p + '/include')
    lib = mlibpath.group(1) if mlibpath else (p + '/lib')
    rv['include_dirs'] += [inc]
    rv['library_dirs'] += [lib]
    return rv

def get_gsl_config_win():
    '''Return dictionary with paths to GSL library, windwows version.
       This version is installed with conda.
    '''
    conda_prefix = os.environ['CONDA_PREFIX']
    inc = os.path.join(conda_prefix, 'Library', 'include')
    lib = os.path.join(conda_prefix, 'Library', 'lib')
    rv = {'include_dirs': [], 'library_dirs': []}
    rv['include_dirs'] += [inc]
    rv['library_dirs'] += [lib]
    return rv

# ----------------------------------------------------------------------------

# compile and link options
define_macros = []
os_name = os.name
if os_name == 'nt':
    gcfg = get_gsl_config_win()
else:
    gcfg = get_gsl_config()
include_dirs = [MYDIR] + gcfg['include_dirs']
library_dirs = []
libraries = []
extra_objects = []
extra_compile_args = []
extra_link_args = []

compiler_type = get_compiler_type()
if compiler_type in ("unix", "cygwin", "mingw32"):
    extra_compile_args = ['-std=c++11', '-Wall', '-Wno-write-strings',
                          '-O3', '-funroll-loops', '-ffast-math']
    extra_objects += ((p + '/libgsl.a') for p in gcfg['library_dirs'])
elif compiler_type == "msvc":
    define_macros += [('_USE_MATH_DEFINES', None)]
    extra_compile_args = ['/EHs']
    libraries += ['gsl']
    library_dirs += gcfg['library_dirs']
# add optimization flags for other compilers if needed


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
        include_dirs = include_dirs,
        libraries = libraries,
        library_dirs = library_dirs,
        define_macros = define_macros,
        extra_compile_args = extra_compile_args,
        extra_link_args = extra_link_args,
        extra_objects = extra_objects,
)


with open(os.path.join(MYDIR, 'README.rst')) as fp:
    long_description = fp.read()

# define distribution
setup_args = dict(
    name = 'diffpy.pdffit2',
    version = versiondata.get('DEFAULT', 'version'),
    packages = find_packages(),
    test_suite = 'diffpy.pdffit2.tests',
    ext_modules = [pdffit2module],
    include_package_data = True,
    install_requires = [
        'six',
        'diffpy.structure>=3',
    ],
    zip_safe = False,

    author = 'Simon J.L. Billinge',
    author_email = 'sb2896@columbia.edu',
    maintainer = 'Pavol Juhas',
    maintainer_email = 'pavol.juhas@gmail.com',
    url = 'https://github.com/diffpy/diffpy.pdffit2',
    description = 'PDFfit2 - real space structure refinement program.',
    long_description = long_description,
    long_description_content_type = 'text/x-rst',
    license = 'BSD',
    keywords = 'PDF structure refinement',
    classifiers = [
        # List of possible values at
        # http://pypi.python.org/pypi?:action=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: C++',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)

if __name__ == '__main__':
    setup(**setup_args)

# End of file
