#CPP= g++-4.5
CPP= g++
CFLAGS=-O -O1 -O2 -O3  -Wl,--no-as-needed
#CFLAGS= -D ITERATION
#DFLAGS= -D DEBUG -g
DFLAGS=-g
PREFIX=/opt
#PREFIX=/home/share

SOURCES= main.o \
	 map.o \
	 simSystem.o \
	 ULAGENT.o \
	 CORRE_MA_OPE.o\
	 clusterStructure.o \
	 mapFactory.o \
	 csFactory.o \
	 utility.o \
	 scheduler.o \
	 MyTMINLP.o \
	 branchBoundScheduler.o \
	 maxSNRScheduler.o \
	 schedFactory.o \
	 simulator.o 
OBJECTS=$(SOURCES:.cpp=.o)

INCLUDE_PATH=-I$(PREFIX)/boost1550/include \
	     -I$(PREFIX)/armadillo40002/include \
	     -I$(PREFIX)/eigen320/include\
	     -I$(PREFIX)/bonmin174/include

LIB_PATH=-L$(PREFIX)/boost1550/lib \
	 -L$(PREFIX)/armadillo40002/lib \
	 -L$(PREFIX)/openblas/lib \
	 -L$(PREFIX)/lapack/lib \
	 -L$(PREFIX)/bonmin174/lib

LIB_FLAGS=-lopenblas \
	  -llapack \
	  -ldl\
	  -larmadillo \
	  -lboost_program_options\
	  -lboost_regex \
	  -lboost_system \
	  -lboost_filesystem\
	  -lboost_math_c99 \
	  -lboost_math_tr1 \
	  -lipopt\
	  -lbonmin\
	  -lbonminampl\
	  -lipoptamplinterface\
	  -lCbc\
	  -lCbcSolver\
	  -lCgl\
	  -lClp\
	  -lCoinUtils\
	  -lOsi\
	  -lOsiCbc\
	  -lOsiClp\
	  -lipopt\
	  -lcoinhsl\
	  -lcoinasl\
	  -lcoinmumps\
	  -lcoinmetis\
	  -lcoinblas\
	  -lcoinlapack\
	  -lz\
	  -lbz2

all: simulator 

simulator: $(OBJECTS)
	@echo Building $@
	@-$(CPP)  $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH)  $^ -o $@ $(LIB_FLAGS)

%.o: %.cpp
	@echo Making $@
	@-$(CPP) -c $(CFLAGS) $(DFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $< -o $@ $(LIB_FLAGS) 
clean:
	rm simulator *.o
