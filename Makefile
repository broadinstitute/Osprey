CXX=		g++
CXXFLAGS=	-g -O3 -fopenmp -Wall
CPPFLAGS=	-I$(HTSLIBDIR)/include/htslib -Iinclude -I$(BOOST_DIR)/include
INCLUDES=	
OBJS=		
EXES=		$(BINDIR)/computeIBS $(BINDIR)/phaseImpMissing $(BINDIR)/osprey
BOOSTLIBS=	$(BOOST_DIR)/lib/libboost_iostreams.a
#LIBS=		-L$(HTSLIBDIR)/lib -lhts -Wl,-Bstatic $(BOOSTLIBS) -lz -Wl,-Bdynamic
#LIBS=		$(HTSLIBDIR)/lib/libhts.a -Wl,-Bstatic $(BOOSTLIBS) -Wl,-Bdynamic -lcurl -lz
LIBS=		$(HTSLIBDIR)/lib/libhts.a ${BOOSTLIBS} -lcurl -lz -llzma -lbz2

BINDIR=		bin
BOOST_DIR=	/humgen/cnp04/bobh/boost_1_58_0/install
HTSLIBDIR=	/humgen/cnp04/bobh/htslib


.SUFFIXES:	.cpp .o
.PHONY:		all clean

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@


all:		$(EXES)

$(BINDIR)/osprey: src/Osprey.o src/OspreyParams.o src/VCFReader.o src/VCFWriter.o src/Variant.o src/PhaseImpMissing.o src/PhaseImpCore.o src/format.o src/timestamp.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

$(BINDIR)/computeIBS: src/computeIBS.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

$(BINDIR)/phaseImpMissing: src/phaseImpMissing.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -f src/*.o a.out $(EXES)
