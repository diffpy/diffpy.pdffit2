########################################################################
# Common targets:     all install install-lib install-scripts
#
# User variables:
#
#   PYTHON_INCLUDE    path to Python include files used for compilation
#
# $Id$
########################################################################

ifndef PYTHON_INCLUDE
PYTHON_INCLUDE = $(shell python -c 'import sys; \
		 print sys.prefix + "/include/python" + sys.version[:3]')
endif

########################################################################

INCLUDE = \
	  -I$(PYTHON_INCLUDE) \
	  -Ilibpdffit2             \
	  -Ipdffit2module          \
	  -Ibuild -I.

DEFINES := $(shell python -c 'import setup_args; setup_args.printDefines()')

GSLLIBS := $(shell gsl-config --libs)

OPTIMFLAGS = -O3 -Wall -Wno-write-strings -funroll-loops -ffast-math -fPIC
DEBUGFLAGS = -gstabs+ -Wall -fPIC

ifdef DEBUG
CPPFLAGS = $(DEBUGFLAGS) $(INCLUDE) $(DEFINES)
else
CPPFLAGS = $(OPTIMFLAGS) $(INCLUDE) $(DEFINES)
endif

all: build build/pdffit2module.so

clean:
	rm -rf -- build

OBJS = \
    build/bindings.o \
    build/misc.o \
    build/pdffit2module.o \
    build/pyexceptions.o \
    build/Atom.o \
    build/LocalPeriodicTable.o \
    build/OutputStreams.o \
    build/PeriodicTable.o \
    build/PointsInSphere.o \
    build/StringUtils.o \
    build/fit.o \
    build/gaussj.o \
    build/metric.o \
    build/nrutil.o \
    build/output.o \
    build/parser.o \
    build/pdf.o \
    build/pdffit.o \
    build/pdflsmin.o \
    build/scatlen.o \
    build/stru.o \


build:
	mkdir $@

build/pdffit2module.so: $(OBJS)
	g++ -o $@ -shared $(OBJS) $(GSLLIBS)

build/%.o : libpdffit2/%.cc
	g++ -c $(CPPFLAGS) -o $@ $<

build/%.o : pdffit2module/%.cc
	g++ -c $(CPPFLAGS) -o $@ $<

# End of file
