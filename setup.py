#!/usr/bin/env python

# Installation script for diffpy.pdffit2

"""PDFfit2 - real space structure refinement engine

Packages:   diffpy.pdffit2
Scripts:    pdffit2
"""

import sys
import os
import warnings
from setuptools import setup, find_packages
from setuptools import Extension

# Use this version when git data are not available, like in git zip archive.
# Update when tagging a new release.
FALLBACK_VERSION = '1.1.post0'

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
MYDIR = os.path.dirname(os.path.abspath(__file__))
versioncfgfile = os.path.join(MYDIR, 'diffpy/pdffit2/version.cfg')
gitarchivecfgfile = versioncfgfile.replace('version.cfg', 'gitarchive.cfg')

def gitinfo():
    from subprocess import Popen, PIPE
    kw = dict(stdout=PIPE, cwd=MYDIR)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %at %ai'], **kw)
    glog = proc.stdout.read()
    rv = {}
    rv['version'] = '.post'.join(desc.strip().split('-')[:2]).lstrip('v')
    rv['commit'], rv['timestamp'], rv['date'] = glog.strip().split(None, 2)
    return rv


def getversioncfg():
    import re
    from ConfigParser import RawConfigParser
    vd0 = dict(version=FALLBACK_VERSION, commit='', date='', timestamp=0)
    # first fetch data from gitarchivecfgfile, ignore if it is unexpanded
    g = vd0.copy()
    cp0 = RawConfigParser(vd0)
    cp0.read(gitarchivecfgfile)
    if '$Format:' not in cp0.get('DEFAULT', 'commit'):
        g = cp0.defaults()
        mx = re.search(r'\btag: v(\d[^,]*)', g.pop('refnames'))
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
        cp.write(open(versioncfgfile, 'w'))
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


def get_gsl_prefix():
    '''Return prefix directory to the gsl library.
    '''
    global _gsl_prefix
    if _gsl_prefix is not None:
        return _gsl_prefix
    gslcfgpaths = [os.path.join(p, 'gsl-config')
            for p in ([MYDIR] + os.environ['PATH'].split(os.pathsep))]
    gslcfgpaths = [p for p in gslcfgpaths if os.path.isfile(p)]
    if not gslcfgpaths:
        warnings.warn("Cannot find gsl-config in MYDIR nor in system PATH.")
        _gsl_prefix = ''
        return get_gsl_prefix()
    gslcfg = gslcfgpaths[0]
    for line in open(gslcfg, 'r'):
        if line.startswith('prefix='):
            _gsl_prefix = line[7:].strip()
            break
    else:
        emsg = "Cannot find 'prefix=' line in {}."
        raise RuntimeError(emsg.format(gslcfg))
    return get_gsl_prefix()
_gsl_prefix = None

# ----------------------------------------------------------------------------

# compile and link options
define_macros = []
include_dirs = [MYDIR, get_gsl_prefix() + '/include']
extra_objects = [get_gsl_prefix() + '/lib/libgsl.a']
extra_compile_args = []
extra_link_args = []

compiler_type = get_compiler_type()
if compiler_type in ("unix", "cygwin", "mingw32"):
    extra_compile_args = ['-Wall', '-Wno-write-strings',
            '-O3', '-funroll-loops', '-ffast-math']
elif compiler_type == "msvc":
    define_macros += [('_USE_MATH_DEFINES', None)]
    extra_compile_args = ['/EHs']
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
        define_macros = define_macros,
        extra_compile_args = extra_compile_args,
        extra_link_args = extra_link_args,
        extra_objects = extra_objects,
)

# define distribution
setup_args = dict(
        name = 'diffpy.pdffit2',
        version = versiondata.get('DEFAULT', 'version'),
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        test_suite = 'diffpy.pdffit2.tests',
        ext_modules = [pdffit2module],
        include_package_data = True,
        install_requires = [
            'diffpy.Structure>=1.2',
        ],
        zip_safe = False,

        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        maintainer = 'Pavol Juhas',
        maintainer_email = 'pavol.juhas@gmail.com',
        url = 'https://github.com/diffpy/diffpy.pdffit2',
        description = 'PDFfit2 - real space structure refinement program.',
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
            'Programming Language :: Python :: 2.5',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

if __name__ == '__main__':
    setup(**setup_args)

# End of file
