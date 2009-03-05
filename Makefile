########################################################################
# This Makefile should be used only for debugging purposes and is
# intended for strong-hearted developers.  End users should use the
# setup.py script instead.
#
# Variables:
#
#   DEBUG             when defined, compile with debug options
#   PYTHON_INCLUDE    path to Python include files used for compilation
#
# $Id$
########################################################################

ifndef PYTHON_INCLUDE
PYTHON_INCLUDE := $(shell python -c \
    'from distutils import sysconfig; print sysconfig.get_python_inc()')
endif

########################################################################

INCLUDE = \
	  -I$(PYTHON_INCLUDE) \
	  -Ilibpdffit2 \
	  -Ipdffit2module \
	  -Ibuild -I. \
	  $(GSL_INCLUDE)

GSL_INCLUDE := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)

COMMONFLAGS = -Wall -Wno-write-strings -fPIC
OPTIMFLAGS = -O3 -funroll-loops -ffast-math $(COMMONFLAGS)
DEBUGFLAGS = -gstabs+ $(COMMONFLAGS)

ifdef DEBUG
CPPFLAGS = $(DEBUGFLAGS) $(INCLUDE) $(DEFINES)
else
CPPFLAGS = $(OPTIMFLAGS) $(INCLUDE) $(DEFINES)
endif

all: build diffpy/pdffit2/pdffit2.so

clean:
	rm -rf -- build diffpy/pdffit2/pdffit2.so

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

diffpy/pdffit2/pdffit2.so: build/pdffit2module.so
	ln -f $< $@

build/pdffit2module.so: $(OBJS)
	$(CXX) -o $@ -shared $(OBJS) $(GSL_LIBS)

build/%.o : libpdffit2/%.cc
	$(CXX) -c $(CPPFLAGS) -o $@ $<

build/%.o : pdffit2module/%.cc
	$(CXX) -c $(CPPFLAGS) -o $@ $<

# End of file
