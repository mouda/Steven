CPP= g++
LIB_PATH= ../commonLibrary/CORRE_MA_OPE.o \
	  ../commonLibrary/SimSystem.o \
	  ../commonLibrary/timesearch.o \
	  ../commonLibrary/TimeStamp.o \
	  ../commonLibrary/ULAGENT.o \
	  ../commonLibrary/ULConstraintSolver.o \
	  ../commonLibrary/ULCS1b.o 
all: DirectAcessPowerMax_MC

DirectAcessPowerMax_MC: 	main.o 
	$(CPP)  $(LIB_PATH) $(CFLAGS) `pkg-config opencv --cflags --libs`  main.o -o DirectAcessPowerMax_MC

main.o:	 main.cpp 
	$(CPP) $(CFLAGS) -c   main.cpp -o main.o
clean:
	rm DirectAcessPowerMax_MC *.o

