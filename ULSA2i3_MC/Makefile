CPP= g++
CFLAGS=-O -O1 -O2 -O3
#DFLAGS= -D DEBUG
LIB_PATH= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSA2i3_MC

ULSA2i3_MC: 	main.o ULSA2i3_MC.o 
	$(CPP)  $(LIB_PATH) $(CFLAGS)  `pkg-config opencv --cflags --libs` $^  -o $@ 

%.o: %.cpp
	$(CPP) -c $(CFLAGS) $(DFLAGS) $< -o $@
clean:
	rm ULSA2i3_MC *.o
