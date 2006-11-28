# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2005  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Notes:
#
# mm STATIC=1
#
# builds a static (.a) instead of default shared library.  Python module is
# then standalone, but needs to be rebuilt after every change in C++ code.
#
# When environment variable TARGET contains debug, optimizations are turned off
# and library is build with extended debugging info (-gstabs+).

include local.def

PROJECT = pdffit2
PACKAGE = libpdffit2

# set environment variable TARGET=debug to compile debugging info
ifneq (,$(findstring debug,$(TARGET)))
    DEV_CXX_FLAGS = -gstabs+
else
    DEV_CXX_FLAGS = -O3 -ffast-math
endif
PROJ_SAR = $(BLD_LIBDIR)/$(PACKAGE).$(EXT_SAR)
PROJ_DLL = $(BLD_BINDIR)/$(PACKAGE).$(EXT_SO)
PROJ_TMPDIR = $(BLD_TMPDIR)/$(PROJECT)/$(PACKAGE)
PROJ_CLEAN += $(PROJ_SAR) $(PROJ_DLL)

PROJ_SRCS = \
    Atom.cc           \
    OutputStreams.cc  \
    PeriodicTable.cc  \
    PointsInSphere.cc \
    StringUtils.cc    \
    fit.cc            \
    gaussj.cc         \
    math.cc           \
    metric.cc         \
    nrutil.cc         \
    output.cc         \
    parser.cc         \
    pdf.cc            \
    pdffit.cc         \
    pdflsmin.cc       \
    scatlen.cc        \
    stru.cc           \


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build the library

ifdef STATIC
## static link of C code into python interface, no need for export
## watch out when building pdffit2module.so -
## it does not check dependence on PROJ_CXX_LIB
PROJ_CXX_LIB = $(BLD_LIBDIR)/$(PACKAGE).a
PROJ_CLEAN += $(PROJ_CXX_LIB)
all: $(PROJ_CXX_LIB)
else
## link C code into shared library, this needs to be exported
all: $(PROJ_SAR) export
endif

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ifeq (Win32, ${findstring Win32, $(PLATFORM_ID)})

# build the shared object
$(PROJ_SAR): product_dirs $(PROJ_OBJS)
	$(CXX) $(LCXXFLAGS) -o $(PROJ_DLL) \
	-Wl,--out-implib=$(PROJ_SAR) $(PROJ_OBJS)

# export
export:: export-headers export-libraries export-binaries

else

# build the shared object
$(PROJ_SAR): product_dirs $(PROJ_OBJS)
	$(CXX) $(LCXXFLAGS) -o $(PROJ_SAR) $(PROJ_OBJS)

# export
export:: export-headers export-libraries

endif

EXPORT_HEADERS = \
    pdffit.h   \


EXPORT_LIBS = $(PROJ_SAR)
EXPORT_BINS = $(PROJ_DLL)


# version
# $Id$

#
# End of file
