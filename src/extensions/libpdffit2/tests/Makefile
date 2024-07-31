########################################################################
# Targets:  alltests unittests
#
# Comments: This Makefile is for building and running C++ unit tests
# on the libpdffit2 sources.  Intended for developers, not for the end
# users who just want to build pdffit2 from sources.
#
########################################################################


# compiler options
CPPFLAGS += -I../ $(shell cppunit-config --cflags)
LDLIBS += $(shell cppunit-config --libs)


# default target
all: unittest


# build unit test driver
alltests: alltests.cc \
    TestPointsInSphere.cc ../PointsInSphere.cc


# execute unit tests
unittest: alltests
	./alltests


# clean up any generated files
clean:
	rm -f alltests
