CPP= g++
#DFLAGS= -D DEBUG
LIB_PATH= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSA_DC_Estimate

ULSA_DC_Estimate: 	main.o ULSA_DC_Estimate.o 
	$(CPP) $(LIB_PATH) $(CFLAGS) $(DFLAGS) `pkg-config opencv --cflags --libs` $^ -o $@

%.o: %.cpp
	$(CPP) -c $(CFLAGS) $(DFLAGS) $< -o $@

clean:
	rm ULSA_DC_Estimate *.o
