bin_PROGRAMS = simulator
INCLUDES =-I../lib
#scheduler_CFLAGS = 
LIBS = -lopenblas \
       -llapack \
       -ldl\
       -larmadillo \
       -lboost_program_options\
       -lboost_regex \
       -lboost_system \
       -lboost_filesystem\
       -lboost_math_c99 \
       -lboost_math_tr1 \
       -lipopt\
       -lbonmin\
       -lbonminampl\
       -lipoptamplinterface\
       -lCbc\
       -lCbcSolver\
       -lCgl\
       -lClp\
       -lCoinUtils\
       -lOsi\
       -lOsiCbc\
       -lOsiClp\
       -lipopt\
       -lcoinhsl\
       -lcoinasl\
       -lcoinmumps\
       -lcoinmetis\
       -lcoinblas\
       -lcoinlapack\
       -lz\
       -lbz2
AM_LDFLAGS = -lopenblas \
	     -llapack \
	     -ldl\
	     -larmadillo \
	     -lboost_program_options\
	     -lboost_regex \
	     -lboost_system \
	     -lboost_filesystem\
	     -lboost_math_c99 \
	     -lboost_math_tr1 \
	     -lipopt\
	     -lbonmin\
	     -lbonminampl\
	     -lipoptamplinterface\
	     -lCbc\
	     -lCbcSolver\
	     -lCgl\
	     -lClp\
	     -lCoinUtils\
	     -lOsi\
	     -lOsiCbc\
	     -lOsiClp\
	     -lipopt\
	     -lcoinhsl\
	     -lcoinasl\
	     -lcoinmumps\
	     -lcoinmetis\
	     -lcoinblas\
	     -lcoinlapack\
	     -lz\
	     -lbz2
simulator_SOURCES = main.cpp \
		    map.h map.cpp \
		    simSystem.h simSystem.cpp \
		    ULAGENT.h ULAGENT.cpp \
		    CORRE_MA_OPE.h CORRE_MA_OPE.cpp\
		    clusterStructure.h clusterStructure.cpp \
		    mapFactory.h mapFactory.cpp \
		    csFactory.h csFactory.cpp \
		    utility.h utility.cpp \
		    scheduler.h scheduler.cpp \
		    MyTMINLP.hpp MyTMINLP.cpp \
		    branchBoundScheduler.h branchBoundScheduler.cpp \
		    maxSNRScheduler.h maxSNRScheduler.cpp \
		    greedyPhysical.h greedyPhysical.cpp \
		    schedFactory.h schedFactory.cpp \
		    slot.h slot.cpp \
		    event.h event.cpp \
		    fileHandler.h fileHandler.cpp\
		    simulator.h simulator.cpp 

if DEBUG
AM_CFLAGS = -g3 -O0
AM_CXXFLAGS = -g3 -O0
else
AM_CFLAGS = -O2
AM_CXXFLAGS = -O2
endif
