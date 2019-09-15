OBJECTS := FragmentBlocks.o ReadBlockProcessor.o CoverageBlock.o ReadBlockProcessor_CoverageBlocks.o ReadBlockProcessor_OutputBAM.o BAM2blocks.o IRFinder.o crc32.o
SOURCES=$(wildcard *.cpp)
LDFLAGS :=
LDFLAGS_static := -static -static-libgcc
LDFLAGS_GDB :=
OPTIMFLAGS :=
OPTIMFLAGS1 :=
# below flags make little difference.
#OPTIMFLAGS=-ffast-math
#OPTIMFLAGS=-fforce-addr -funsafe-loop-optimizations -ftree-loop-linear -ftree-vectorize 
#OPTIMFLAGS=-fforce-addr -funsafe-loop-optimizations -ftree-vectorize 
# unroll-loops slows the program ~6%.
#OPTIMFLAGS1=-funroll-loops  -fprefetch-loop-arrays
COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(shell echo `date` `hostname --fqdn`:`pwd`)"'
CCFLAGS_MAIN := -pipe -std=c++0x -O3    -Wall -Wextra -fopenmp $(COMPTIMEPLACE) $(OPTIMFLAGS) $(OPTIMFLAGS1)
CCFLAGS_GDB := -pipe  -std=c++0x -g -O0 -Wall -Wextra -fopenmp $(COMPTIMEPLACE)
CC :=g++
GCC:=gcc

%.o : %.cpp
	$(CC) -c $(CCFLAGS) $<

all: IRFinderStatic

clean :
	rm -f *.o irfinder Depend.list

ifneq ($(MAKECMDGOALS),clean)
Depend.list: $(SOURCES)
	/bin/rm -f ./Depend.list
	$(CC) $(CCFLAGS_MAIN) -MM $^ >> Depend.list
include Depend.list
endif

IRFinder : CCFLAGS=$(CCFLAGS_MAIN)
IRFinder : $(OBJECTS)
	$(CC) -o irfinder $(CCFLAGS) $(OBJECTS) $(LDFLAGS)

IRFinderStatic : CCFLAGS=$(CCFLAGS_MAIN)
IRFinderStatic : $(OBJECTS)
	$(CC) -o irfinder $(CCFLAGS) $(OBJECTS) $(LDFLAGS_static)

gdb : CCFLAGS= $(CCFLAGS_GDB)
gdb : $(OBJECTS)
	$(CC) -o irfinder $(CCFLAGS_GDB) $(OBJECTS) $(LDFLAGS_GDB) 


