#CPP= g++-4.5
CPP= g++
#CFLAGS=-O -O1 -O2 -O3
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
	 MyNLP.o \
	 branchBoundScheduler.o \
	 maxSNRScheduler.o \
	 schedFactory.o \
	 simulator.o 
OBJECTS=$(SOURCES:.cpp=.o)

INCLUDE_PATH=-I$(PREFIX)/boost1550/include \
	     -I$(PREFIX)/armadillo39202/include \
	     -I$(PREFIX)/eigen320/include\
	     -I$(PREFIX)/CoinIpopt3117/include

LIB_PATH=-L$(PREFIX)/boost1550/lib \
	 -L$(PREFIX)/armadillo39202/lib \
	 -L$(PREFIX)/openblas/lib \
	 -L$(PREFIX)/lpack342/lib \
	 -L$(PREFIX)/CoinIpopt3117/lib

LIB_FLAGS=-larmadillo \
	  -lopenblas \
	  -llapack \
	  -lboost_regex \
	  -lboost_system \
	  -lboost_filesystem\
	  -lboost_filesystem-mt\
	  -lboost_math_c99 \
	  -lboost_math_tr1 \
	  -lipopt\
	  -ldl\
	  -lcoinhsl\
	  -lcoinmumps\
	  -lcoinmetis\
	  -lcoinblas\
	  -lcoinlapack\

all: simulator 

simulator: $(OBJECTS)
	@echo Building $@
	@-$(CPP)  $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH)  $^ -o $@ $(LIB_FLAGS)

%.o: %.cpp
	@echo Making $@
	@-$(CPP) -c $(CFLAGS) $(DFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $< -o $@ $(LIB_FLAGS) 
clean:
	rm simulator *.o
