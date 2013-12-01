CPP= g++
CFLAGS=-O -O1 -O2 -O3
#CFLAGS= -D ITERATION
#DFLAGS= -D DEBUG -g
PREFIX=/opt
#PREFIX=/home/share

SOURCES=main.o ULSA4b8_DC.o 
OBJECTS=$(SOURCES:.cpp=.o)

INCLUDE_PATH=-I$(PREFIX)/boost1550/include \
	     -I$(PREFIX)/armadillo39202/include \
	     -I$(PREFIX)/eigen320/include

LIB_PATH=-L$(PREFIX)/boost1550/lib \
	 -L$(PREFIX)/armadillo39202/lib \
	 -L$(PREFIX)/openblas/lib \
	 -L$(PREFIX)/lpack342/lib 
LIB_FLAGS=-lboost_math_c99 \
	  -lboost_math_tr1 \
	  -larmadillo \
	  -lopenblas \
	  -llapack  

COMM_UTILITY= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSA4b8_DC

ULSA4b8_DC: $(OBJECTS)
	@echo Building $@
	@-$(CPP)  $(COMM_UTILITY) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) $^  -o $@ 

%.o: %.cpp
	@echo Making $@
	@-$(CPP) -c $(CFLAGS) $(DFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) $< -o $@
clean:
	rm ULSA4b8_DC *.o
