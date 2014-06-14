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
		    csPowerUpdater.h csPowerUpdater.cpp\
		    ULSA4b7_DC.h ULSA4b7_DC.cpp\
		    tier1NLP.h tier1NLP.cpp\
		    minPowerSACluster.h minPowerSACluster.cpp\
		    minPowerCsFactory.h minPowerCsFactory.cpp\
		    kmeansCsFactory.h kmeansCsFactory.cpp\
		    minResCsFactory.h minResCsFactory.cpp\
		    scheduler.h scheduler.cpp \
		    MyTMINLP.hpp MyTMINLP.cpp \
		    branchBoundScheduler.h branchBoundScheduler.cpp \
		    maxSNRScheduler.h maxSNRScheduler.cpp \
		    greedyPhysical.h greedyPhysical.cpp \
		    bruteForceScheduler.h bruteForceScheduler.cpp\
		    nonSimplifiedScheduler.h nonSimplifiedScheduler.cpp\
		    minPowerScheduler.h minPowerScheduler.cpp\
		    minPowerMILP.h minPowerMILP.cpp\
		    SASolver.h SASolver.cpp\
		    powerUpdater.h powerUpdater.cpp\
		    schedFactory.h schedFactory.cpp \
		    slot.h slot.cpp \
		    simulator.h simulator.cpp 

LOCAL_INCLUDE = -I../lib/
if DEBUG
AM_CFLAGS = $(BONMIN_CPPFLAGS) $(BONMIN_LIBS)  $(LOCAL_INCLUDE)  -I../ -g3 -O0
AM_CXXFLAGS = $(BONMIN_CPPFLAGS) $(BONMIN_LIBS) -I../ -g3 -O0
else
AM_CFLAGS = $(BONMIN_CPPFLAGS) $(BONMIN_LIBS)  -I../ -O2
AM_CXXFLAGS = $(BONMIN_CPPFLAGS) $(BONMIN_LIBS)  -I../ -O2
endif
simulator_LDADD = $(BONMIN_CPPFLAGS) $(BONMIN_LIBS) ../lib/libfileHandler.la \
		  ../lib/libfileHandler.la \
                  ../lib/libutility.la \
                  ../lib/libmapFactory.la \
                  ../lib/libcsFactory.la \
                  ../lib/libCORRE_MA_OPE.la \
                  ../lib/libULAGENT.la \
                  ../lib/libsimSystem.la \
		  ../lib/libtimesearch.la\
		  ../lib/libTimeStamp.la \
		  ../lib/libULConstraintSolver.la\
		  ../lib/libULCS1b.la
