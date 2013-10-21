CPP=g++
all: CORRE_MA_OPE.o SimSystem.o timesearch.o TimeStamp.o ULAGENT.o ULConstraintSolver.o ULCS1b.o

CORRE_MA_OPE.o:		CORRE_MA_OPE.cpp CORRE_MA_OPE.h 
	$(CPP) $(CFLAGS) -c $< -o $@ `pkg-config opencv --cflags --libs`

SimSystem.o: 		SimSystem.cpp SimSystem.h 
	$(CPP) $(CFLAGS) -c $< -o $@

timesearch.o: 		timesearch.cpp timesearch.h 
	$(CPP) $(CFLAGS) -c $< -o $@

TimeStamp.o: 		TimeStamp.cpp TimeStamp.h 
	$(CPP) $(CFLAGS) -c $< -o $@

ULAGENT.o: 		ULAGENT.cpp ULAGENT.h 
	$(CPP) $(CFLAGS) -c $< -o $@

ULCS1b.o:		ULCS1b.cpp ULCS1b.h 
	$(CPP) $(CFLAGS) -c $< -o $@

ULConstraintSolver.o: 	ULConstraintSolver.cpp ULConstraintSolver.h  
	$(CPP) $(CFLAGS) -c $< -o $@

clean:
	rm *.o 
