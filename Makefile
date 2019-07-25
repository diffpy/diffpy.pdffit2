########################################################################
# This Makefile should be used only for debugging purposes and is
# intended for strong-hearted developers.  End users should use the
# setup.py script instead.
#
# Variables:
#
#   DEBUG             when defined, compile with debug options
#   USESHARED         when defined, link with shared gsl library
#   PYTHON_INCLUDE    path to Python include files used for compilation
#
########################################################################

ifndef PYTHON_INCLUDE
PYTHON_INCLUDE := $(shell python -c \
    'from distutils import sysconfig; print(sysconfig.get_python_inc())')
endif

########################################################################

INCLUDE = \
	  -I$(PYTHON_INCLUDE) \
	  -Ilibpdffit2 \
	  -Ipdffit2module \
	  -Ibuild -I. \
	  $(GSL_INCLUDE)

GSL_INCLUDE = $(shell gsl-config --cflags)

GSL_LIBS_CONFIG = $(shell gsl-config --libs-without-cblas)
ifdef USESHARED
GSL_LIBS = $(GSL_LIBS_CONFIG)
else
GSL_LIBS = $(subst -L,,$(filter -L%, $(GSL_LIBS_CONFIG)))/libgsl.a
endif

ifneq (,$(findstring darwin,$(OSTYPE)))
LDFLAGS = -bundle -undefined dynamic_lookup
else
LDFLAGS = -shared
endif

COMMONFLAGS = -std=c++11 -Wall -Wno-write-strings -fPIC
OPTIMFLAGS = -O3 -funroll-loops -ffast-math $(COMMONFLAGS)
DEBUGFLAGS = -g $(COMMONFLAGS)

ifdef DEBUG
CPPFLAGS = $(DEBUGFLAGS) $(INCLUDE) $(DEFINES)
else
CPPFLAGS = $(OPTIMFLAGS) $(INCLUDE) $(DEFINES)
endif

PYTHON_XY := $(shell python -c 'import sys; print(sys.version[:3])')
BDIR := build/$(PYTHON_XY)

.PHONY: all module clean test

all: module

module: $(BDIR)/pdffit2module.so
	@test $< -ef diffpy/pdffit2/pdffit2.so || \
	    ln -vf $< diffpy/pdffit2/pdffit2.so

clean:
	rm -rf -- build/$(PYTHON_XY) diffpy/pdffit2/pdffit2.so

test: module
	python -m diffpy.pdffit2.tests.run

OBJS = \
    $(BDIR)/bindings.o \
    $(BDIR)/misc.o \
    $(BDIR)/pdffit2module.o \
    $(BDIR)/pyexceptions.o \
    $(BDIR)/Atom.o \
    $(BDIR)/LocalPeriodicTable.o \
    $(BDIR)/OutputStreams.o \
    $(BDIR)/PeriodicTable.o \
    $(BDIR)/PointsInSphere.o \
    $(BDIR)/StringUtils.o \
    $(BDIR)/fit.o \
    $(BDIR)/gaussj.o \
    $(BDIR)/metric.o \
    $(BDIR)/nrutil.o \
    $(BDIR)/output.o \
    $(BDIR)/parser.o \
    $(BDIR)/pdf.o \
    $(BDIR)/pdffit.o \
    $(BDIR)/pdflsmin.o \
    $(BDIR)/scatlen.o \
    $(BDIR)/stru.o \


$(BDIR):
	mkdir -p $@

$(BDIR)/pdffit2module.so: $(BDIR) $(OBJS)
	$(CXX) -o $@ $(LDFLAGS) $(OBJS) $(GSL_LIBS)

$(BDIR)/%.o : libpdffit2/%.cc
	$(CXX) -c $(CPPFLAGS) -o $@ $<

$(BDIR)/%.o : pdffit2module/%.cc
	$(CXX) -c $(CPPFLAGS) -o $@ $<

# End of file
