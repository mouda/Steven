CPP= g++
#DFLAGS= -D DEBUG
LIB_PATH= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: ULSA4b3_DC

ULSA4b3_DC: 	main.o ULSA4b3_DC.o 
	$(CPP)  $(LIB_PATH) $(CFLAGS) $(DFLAGS) `pkg-config opencv --cflags --libs`  main.o ULSA4b3_DC.o -o ULSA4b3_DC

main.o: 	../commonLibrary/ULSAOutputToolSet.h main.cpp 
	$(CPP) $(CFLAGS) $(DFLAGS) -c   main.cpp -o main.o

ULSA4b3_DC.o:  ../commonLibrary/ULAGENT.h ULSA4b3_DC.h ULSA4b3_DC.cpp
	$(CPP) $(CFLAGS) $(DFLAGS) -c  ULSA4b3_DC.cpp -o ULSA4b3_DC.o

clean:
	rm ULSA4b3_DC *.o
