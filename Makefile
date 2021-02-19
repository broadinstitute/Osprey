CXX=		g++
CXXFLAGS=	-g -O3 -fopenmp -Wall
CPPFLAGS=	-Iinclude -I$(BOOST_DIR)/include
INCLUDES=
OBJS=		
EXES=		$(BINDIR)/computeIBS $(BINDIR)/phaseImpMissing
LIBS=		-Wl,-Bstatic -lboost_iostreams -lz -Wl,-Bdynamic -L$(BOOST_DIR)/lib

BINDIR=		bin
BOOST_DIR=	/humgen/cnp04/bobh/boost_1_58_0/install

#
#g++ -g -fopenmp -O3 -Wall \
#    phaseImpMissing.cpp \
#    -Iinclude -I${boostDir}/include \
#    -Wl,-Bstatic -lboost_iostreams -lz -Wl,-Bdynamic \
#    -o phaseImpMissing \
#    -L${boostDir}/lib \
#    || exit 1

.SUFFIXES:	.cpp .o
.PHONY:		all clean

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:		$(EXES)

$(BINDIR)/computeIBS: src/computeIBS.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

$(BINDIR)/phaseImpMissing: src/phaseImpMissing.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -f src/*.o a.out $(EXES)
