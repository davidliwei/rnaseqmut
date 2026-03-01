CXXFLAGS += -g -O2
CC = gcc
CPP = g++

BAMTOOLS_BUILD_DIR = bamtools/build
BAMTOOLS_LIB = $(BAMTOOLS_BUILD_DIR)/src/api/libbamtools.a

INCLUDES = -I./ -Ibamtools/src -Ibamtools/src/api -Ibamtools/src/utils -I$(BAMTOOLS_BUILD_DIR)/src/api
LIBRARY = -lz
PROGRAM = ../bin/rnaseqmut

.PHONY: all clean

all: $(PROGRAM)

$(PROGRAM): bamcallmut.o bamalignfunc.o parseargs.o mutmap.o givenrefmut.o refio.o $(BAMTOOLS_LIB)
	$(CPP) $(CXXFLAGS) $(INCLUDES) $^ $(LIBRARY) -o $@

$(BAMTOOLS_LIB):
	cmake -S bamtools -B $(BAMTOOLS_BUILD_DIR)
	cmake --build $(BAMTOOLS_BUILD_DIR) --target BamTools

.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -rf *.o $(PROGRAM)
