CXXFLAGS+=	-g -O2
CC=		gcc
CPP=		g++
#INCLUDES=	-I./ -Ibamtools -Ibamtools/include -Lbamtools/lib
INCLUDES=	-I./  -Ibamtools/build/src/api/include -Ibamtools/build/src/include -Lbamtools/lib
# LIBRARY=	-lbamtools -lz
#LIBRARY=	bamtools/lib/libbamtools.a -lz
LIBRARY= 	bamtools/build/src/api/libbamtools.a -lz
PROGRAM=	../bin/rnaseqmut

.PHONY:		all clean

all:		${PROGRAM}

${PROGRAM}:	bamcallmut.o bamalignfunc.o  parseargs.o mutmap.o givenrefmut.o refio.o
	${CPP} ${CXXFLAGS}  ${INCLUDES}  $^ ${LIBRARY} -o $@


.SUFFIXES:	.cpp .o

.cpp.o:
	${CC} -c ${CXXFLAGS} ${INCLUDES} $< -o $@

clean:
	rm -rf *.o ${PROGRAM}
