CPP= g++
#DFLAGS= -D DEBUG
LIB_PATH= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSA4b2_DC

ULSA4b2_DC: 	main.o ULSA4b2_DC.o 
	$(CPP)  $(LIB_PATH) $(CFLAGS)  `pkg-config opencv --cflags --libs` $^  -o $@ 

%.o: %.cpp
	$(CPP) -c $(CFLAGS) $(DFLAGS) $< -o $@
clean:
	rm ULSA4b2_DC *.o
