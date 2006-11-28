########################################################################
# Common targets:     all install install-lib install-scripts
#
# User variables:
#
#   PYTHON_VERSION    version of Python used for compilation, e.g. "2.3"
#
#   PYTHON_LIB_PATH   install directory for diffpy.pdffit2 module, must
#                     be add to PYTHONPATH if changed from default value
#
#   BINDIR            install directory for command-line pdffit2
#
# $Id$
########################################################################

ifndef PYTHON_VERSION
PYTHON_VERSION = $(shell python -c 'import sys; print sys.version[:3]')
endif
PYTHON_LIB_PATH = /usr/lib/python$(PYTHON_VERSION)/site-packages
BINDIR = /usr/local/bin

########################################################################

PYTHON = python$(PYTHON_VERSION)
BDIFFPY = build/diffpy
IDIFFPY = $(PYTHON_LIB_PATH)/diffpy

INCLUDE = \
	  -I/usr/include/$(PYTHON) \
	  -Ilibpdffit2             \
	  -Ipdffit2module          \
	  -Ibuild -I.

DEFINES := $(shell $(PYTHON) -c 'import setup; setup.printDefines()')

GSLLIBS := $(shell gsl-config --libs)

OPTIMFLAGS = -O3 -Wall -funroll-loops -ffast-math
DEBUGFLAGS = -gstabs+ -Wall

ifdef DEBUG
CPPFLAGS = $(DEBUGFLAGS) $(INCLUDE) $(DEFINES)
else
CPPFLAGS = $(OPTIMFLAGS) $(INCLUDE) $(DEFINES)
endif
	
OBJS = \
    build/bindings.o \
    build/misc.o \
    build/pdffit2module.o \
    build/pyexceptions.o \
    build/Atom.o \
    build/OutputStreams.o \
    build/PeriodicTable.o \
    build/PointsInSphere.o \
    build/StringUtils.o \
    build/fit.o \
    build/gaussj.o \
    build/math.o \
    build/metric.o \
    build/nrutil.o \
    build/output.o \
    build/parser.o \
    build/pdf.o \
    build/pdffit.o \
    build/pdflsmin.o \
    build/scatlen.o \
    build/stru.o \

PYMODULES = \
    $(BDIFFPY)/pdffit2/__init__.py \
    $(BDIFFPY)/pdffit2/version.py  \
    $(BDIFFPY)/pdffit2/PdfFit.py

all: $(BDIFFPY)/pdffit2 $(BDIFFPY)/pdffit2/pdffit2module.so $(PYMODULES)

clean:
	rm -rf -- build

$(BDIFFPY)/pdffit2/pdffit2module.so: $(OBJS)
	g++ -o $@ -shared $(OBJS) $(GSLLIBS)

$(BDIFFPY)/pdffit2:
	mkdir -p $@
	if [ ! -f $(BDIFFPY)/__init__.py ]; then \
	    $(PYTHON) -c 'import setup; print setup.diffpy__init__code' \
		> $(BDIFFPY)/__init__.py; \
	fi

prepare-diffpy:
	mkdir -p $(IDIFFPY)/pdffit2
	if [ ! -f $(IDIFFPY)/__init__.py ]; then \
	    install -m 644 $(BDIFFPY)/__init__.py $(IDIFFPY) \
	fi


install-lib: prepare-diffpy
	install -m 755 $(BDIFFPY)/pdffit2/pdffit2module.so $(IDIFFPY)/pdffit2
	install -m 644 $(BDIFFPY)/pdffit2/*.py $(IDIFFPY)/pdffit2
	$(PYTHON) \
	    /usr/lib/python$(PYTHON_VERSION)/compileall.py \
	    $(PYTHON_LIB_PATH)/diffpy

install-scripts:
	install -D -m 755 applications/pdffit2 $(BINDIR)/pdffit2

install: install-lib install-scripts

build/%.o : libpdffit2/%.cc
	g++ -c $(CPPFLAGS) -o $@ $<

build/%.o : pdffit2module/%.cc
	g++ -c $(CPPFLAGS) -o $@ $<

build/diffpy/pdffit2/%.py : pdffit2/%.py
	cp -pv -- $< $@

# End of file
