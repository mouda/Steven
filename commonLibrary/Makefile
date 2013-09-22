
CPP=g++
all: CORRE_MA_OPE.a SimSystem.a timesearch.a TimeStamp.a ULAGENT.a ULConstraintSolver.a ULCS1b.a 


CORRE_MA_OPE.a: CORRE_MA_OPE.o
	ar rvs CORRE_MA_OPE.a CORRE_MA_OPE.o
SimSystem.a: SimSystem.o
	ar rvs SimSystem.a SimSystem.o
timesearch.a: timesearch.o
	ar rvs timesearch.a timesearch.o
TimeStamp.a: TimeStamp.o
	ar rvs TimeStamp.a TimeStamp.o
ULAGENT.a: ULAGENT.o
	ar rvs ULAGENT.a ULAGENT.o
ULConstraintSolver.a: ULConstraintSolver.o
	ar rvs ULConstraintSolver.a ULConstraintSolver.o
ULCS1b.a: ULCS1b.o
	ar rvs ULCS1b.a	ULCS1b.o
#ULSAOutputToolSet.a: ULSAOutputToolSet.o
#	ar rvs ULSAOutputToolSet.a ULSAOutputToolSet.o

CORRE_MA_OPE.o:		CORRE_MA_OPE.h CORRE_MA_OPE.cpp
	$(CPP) $(CFLAGS) -c CORRE_MA_OPE.cpp -o CORRE_MA_OPE.o `pkg-config opencv --cflags --libs`

SimSystem.o: 		SimSystem.h SimSystem.cpp
	$(CPP) $(CFLAGS) -c SimSystem.cpp -o SimSystem.o

timesearch.o: 		timesearch.h timesearch.cpp
	$(CPP) $(CFLAGS) -c timesearch.cpp -o timesearch.o

TimeStamp.o: 		TimeStamp.h TimeStamp.cpp
	$(CPP) $(CFLAGS) -c TimeStamp.cpp -o TimeStamp.o

ULAGENT.o: 		ULAGENT.h ULAGENT.cpp
	$(CPP) $(CFLAGS) -c ULAGENT.cpp -o ULAGENT.o

ULCS1b.o:		ULCS1b.h ULCS1b.cpp
	$(CPP) $(CFLAGS) -c ULCS1b.cpp -o ULCS1b.o

ULConstraintSolver.o: 	ULConstraintSolver.h ULConstraintSolver.cpp 
	$(CPP) $(CFLAGS) -c ULConstraintSolver.cpp -o ULConstraintSolver.o


#ULSAOutputToolSet.o:	ULSAOutputToolSet.h ULSAOutputToolSet.cpp
#	$(CPP) -c ULSAOutputToolSet.cpp -o ULSAOutputToolSet.o

clean:
	rm *.o *.a
