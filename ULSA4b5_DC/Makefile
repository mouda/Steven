CPP= g++
#CFLAGS=-O -O1 -O2 -O3
CFLAGS= -g -D DEBUG 
#DFLAGS= -D DEBUG

LIB_PATH= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSA4b5_DC

ULSA4b5_DC: 	main.o ULSA4b5_DC.o 
	$(CPP)  $(LIB_PATH) $(CFLAGS)  `pkg-config opencv --cflags --libs` $^  -o $@ 

%.o: %.cpp
	$(CPP) `pkg-config opencv --cflags --libs` -c $(CFLAGS) $(DFLAGS) $< -o $@
clean:
	rm ULSA4b5_DC *.o
