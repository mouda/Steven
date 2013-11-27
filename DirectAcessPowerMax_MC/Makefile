CPP= g++

#PREFIX=/home/share
PREFIX=/opt
INCLUDE_PATH=-I$(PREFIX)/boost1550/include \
	     -I$(PREFIX)/armadillo39202/include
LIB_PATH=-L$(PREFIX)/boost1550/lib \
	 -L$(PREFIX)/armadillo39202/lib \
	 -L$(PREFIX)/openblas/lib \
	 -L$(PREFIX)/lpack342/lib 
LIB_FLAGS=-lboost_math_c99 \
	  -lboost_math_tr1 \
	  -larmadillo 

COMM_UTILITY= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 

all: DirectAcessPowerMax_MC

DirectAcessPowerMax_MC: 	main.o 
	@echo Building $@
	@-$(CPP)  $(COMM_UTILITY) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS)  $< -o $@ 

main.o:	 main.cpp 
	@echo Making $@
	@-$(CPP) $(CFLAGS) -c $< -o $@
clean:
	rm DirectAcessPowerMax_MC *.o

