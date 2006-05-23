# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2004  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = pdffit2
PACKAGE = pdffit2module
MODULE = pdffit2

include std-pythonmodule.def
include local.def

PROJ_CXX_SRCLIB = -lpdffit2

PROJ_SRCS = \
    bindings.cc \
    exceptions.cc \
    misc.cc
    

# version
# $Id: Make.mm,v 1.3 2005/05/18 21:07:20 farrowch Exp $

# End of file
