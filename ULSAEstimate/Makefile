CPP= g++
CFLAGS=-O -O1 -O2 -O3
#CFLAGS= -D ITERATION
#DFLAGS= -D DEBUG -g
#PREFIX=/opt
PREFIX=/home/share
INCLUDE= -I$(PREFIX)/armadillo39202/include/
LIB_PATH = -L$(PREFIX)/armadillo39202/lib/ \
	   -L$(PREFIX)/openblas/lib \
	   -L$(PREFIX)/lpack342/lib/ 
LIB_FLAGS= -larmadillo -lopenblas -llapack  
COMM_UTILITY= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSAEstimate

ULSAEstimate: 	main.o ULSAEstimate.o 
	$(CPP)  $(COMM_UTILITY) $(CFLAGS) $(INCLUDE) $(LIB_PATH) $(LIB_FLAGS) `pkg-config opencv --cflags --libs` $^  -o $@ 

%.o: %.cpp
	$(CPP) -c $(CFLAGS) $(DFLAGS) $(INCLUDE) $(LIB_PATH) $(LIB_FLAGS) $< -o $@
clean:
	rm ULSAEstimate *.o
