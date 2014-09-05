FC=gfortran
CC=g++
CXX=g++
FFLAGS=
CCFLAGS=-fPIC -g -Wall
EXEC=root2lund

ROOTCONFIG   := ${ROOTSYS}/bin/root-config

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)

all: $(EXEC)

root2lund : roottolund.o
	$(CC) $^ $(ROOTCFLAGS) $(ROOTLIBS)  -o $@ 
	$(CC) $^ $(ROOTCFLAGS) $(ROOTLIBS) -DUSE_SOLLUND -o root2sollund

roottolund.o : roottolund.C
	$(CC)  -c $(CCFLAGS) $(ROOTCFLAGS) $< -o $@



#####################################

clean :
	-rm -f *.o $(EXEC)
