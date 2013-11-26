CPP=g++
CFLAGS = -O -O1 -O2 -O3
#CFLAGS = -g

#PREFIX=/home/share
PREFIX=/opt
INCLUDE_PATH=-I$(PREFIX)/boost1550/include \
	     -I$(PREFIX)/armadillo39202/include \
	     -I$(PREFIX)/eigen320/include
LIB_PATH=-L$(PREFIX)/boost1550/lib \
	 -L$(PREFIX)/armadillo39202/lib \
	 -L$(PREFIX)/openblas/lib \
	 -L$(PREFIX)/lpack342/lib 
LIB_FLAGS=-lboost_math_c99 \
	  -lboost_math_tr1 \
	  -larmadillo

all: CORRE_MA_OPE.o SimSystem.o timesearch.o TimeStamp.o ULAGENT.o ULConstraintSolver.o ULCS1b.o

CORRE_MA_OPE.o:		CORRE_MA_OPE.cpp CORRE_MA_OPE.h 
	#$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@ `pkg-config opencv --cflags --libs`
	# opencv lib will be dropped
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@ 

SimSystem.o: 		SimSystem.cpp SimSystem.h 
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@

timesearch.o: 		timesearch.cpp timesearch.h 
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@

TimeStamp.o: 		TimeStamp.cpp TimeStamp.h 
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@

ULAGENT.o: 		ULAGENT.cpp ULAGENT.h 
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@

ULCS1b.o:		ULCS1b.cpp ULCS1b.h 
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@

ULConstraintSolver.o: 	ULConstraintSolver.cpp ULConstraintSolver.h  
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@

clean:
	rm *.o 
