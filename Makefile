INCLUDE = \
	  -I/usr/include/python$(PYTHON_VERSION) \
	  -Ilibpdffit2             \
	  -Ipdffit2module          \
	  -Ibuild -I.

CPPFLAGS = -Wall -O3 \
	   -funroll-loops -fstrict-aliasing -fpic \
	   -ffast-math \
	   $(INCLUDE)
	
OBJS = \
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
    build/lazy.o \
    build/PointsInSphere.o \
    build/bindings.o \
    build/exceptions.o \
    build/misc.o \
    build/pdffit2module.o

all: build build/pdffit2module.so

clean:
	rm -f -- $(OBJS) build/portinfo build/pdffit2module.so
	rmdir build

build/pdffit2module.so: $(OBJS)
	g++ -o $@ -shared -lg2c $(OBJS)

build:
	mkdir build
	touch build/portinfo

build/fit.o: libpdffit2/fit.cc
build/gaussj.o: libpdffit2/gaussj.cc
build/math.o: libpdffit2/math.cc
build/metric.o: libpdffit2/metric.cc
build/nrutil.o: libpdffit2/nrutil.cc
build/output.o: libpdffit2/output.cc
build/parser.o: libpdffit2/parser.cc
build/pdf.o: libpdffit2/pdf.cc
build/pdffit.o: libpdffit2/pdffit.cc
build/pdflsmin.o: libpdffit2/pdflsmin.cc
build/scatlen.o: libpdffit2/scatlen.cc
build/stru.o: libpdffit2/stru.cc
build/lazy.o: libpdffit2/lazy.f
build/PointsInSphere.o: libpdffit2/PointsInSphere.cc
build/bindings.o: pdffit2module/bindings.cc
build/exceptions.o: pdffit2module/exceptions.cc
build/misc.o: pdffit2module/misc.cc
build/pdffit2module.o: pdffit2module/pdffit2module.cc

build/%.o : libpdffit2/%.cc
	g++ -c $(CPPFLAGS) -o $@ $<

build/%.o : libpdffit2/%.f
	g77 -c $(CPPFLAGS) -o $@ $<

build/%.o : pdffit2module/%.cc
	g++ -c $(CPPFLAGS) -o $@ $<
