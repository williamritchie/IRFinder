OBJECTS := TrimReads.o sequenceTools.o trim.o 
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
SVNDEF := -D'SVN_VERSION_COMPILED="magictrim0.1"'
COMPTIMEPLACE := -D'COMPILATION_TIME_PLACE="$(shell echo `date` `hostname --fqdn`:`pwd`)"'
CCFLAGS_MAIN := -pipe -std=c++0x -O3    -Wall -Wextra -fopenmp $(SVNDEF) $(COMPTIMEPLACE) $(OPTIMFLAGS) $(OPTIMFLAGS1)
CCFLAGS_GDB := -pipe  -std=c++0x -g -O0 -Wall -Wextra -fopenmp $(SVNDEF) $(COMPTIMEPLACE)
CC :=g++
GCC:=gcc

%.o : %.cpp
	$(CC) -c $(CCFLAGS) $<

all: trimstatic

clean :
	rm -f *.o trim Depend.list

ifneq ($(MAKECMDGOALS),clean)
Depend.list: $(SOURCES)
	/bin/rm -f ./Depend.list
	$(CC) $(CCFLAGS_MAIN) -MM $^ >> Depend.list
include Depend.list
endif

trim : CCFLAGS=$(CCFLAGS_MAIN)
trim : $(OBJECTS)
	$(CC) -o trim $(CCFLAGS) $(OBJECTS) $(LDFLAGS)

trimstatic : CCFLAGS=$(CCFLAGS_MAIN)
trimstatic : $(OBJECTS)
	$(CC) -o trim $(CCFLAGS) $(OBJECTS) $(LDFLAGS_static)

gdb : CCFLAGS= $(CCFLAGS_GDB)
gdb : $(OBJECTS)
	$(CC) -o trim $(CCFLAGS_GDB) $(OBJECTS) $(LDFLAGS_GDB) 


