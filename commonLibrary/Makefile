PREFIX=/home/share
#PREFIX=/opt

CPP=@-g++
CFLAGS = -O -O1 -O2 -O3 -Wall
#CFLAGS = -g

SOURCES=CORRE_MA_OPE.cpp \
	SimSystem.cpp \
	timesearch.cpp \
	TimeStamp.cpp \
	ULAGENT.cpp \
	ULConstraintSolver.cpp \
	ULCS1b.cpp

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
	  -lboost_filesystem \
	  -lboost_system \
	  -larmadillo

all: $(OBJECTS)

%.o: %.cpp
	@echo Making $@
	$(CPP) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIB_FLAGS) -c $< -o $@ 
clean:
	rm *.o 
