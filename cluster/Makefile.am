bin_PROGRAMS = genCluster 
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
genCluster_SOURCES = main.cpp \
		    myTier1NLP.h myTier1NLP.cpp\
		    ULSA4b7_DC.h ULSA4b7_DC.cpp\
		    minResCsFactory.h minResCsFactory.cpp\
		    kmeansCsFactory.h kmeansCsFactory.cpp\
		    minPowerSACluster.h minPowerSACluster.cpp\
		    minPowerCsFactory.h minPowerCsFactory.cpp\
		    minPowerImageCluster.h minPowerImageCluster.cpp\
		    minPowerImageCsFactory.h minPowerImageCsFactory.cpp\
		    baselineImageCluster.h baselineImageCluster.cpp\
		    baselineImageCsFactory.h baselineImageCsFactory.cpp\
		    nonGuidedSACluster.cpp nonGuidedSACluster.h\
		    nonGuidedCsFactory.cpp nonGuidedCsFactory.h\
		    simulator.h simulator.cpp 

LOCAL_INCLUDE = -I../lib/
if DEBUG
AM_CFLAGS = $(IPOPT_CPPFLAGS) $(LOCAL_INCLUDE)  -I../ -g3 -ggdb -O0
AM_CXXFLAGS = $(IPOPT_CPPFLAGS) $(LOCAL_INCLUDE) -I../ -g3 -ggdb -O0
genCluster_LDFLAGS= $(IPOPT_LIBS) -I../ -g3 -ggdb -O0
else
AM_CFLAGS = $(IPOPT_CPPFLAGS) $(LOCAL_INCLUDE)  -I../ -O2
AM_CXXFLAGS = $(IPOPT_CPPFLAGS) $(LOCAL_INCLUDE)  -I../ -O2
genCluster_LDFLAGS= $(IPOPT_LIBS) -I../ -O2 
endif
genCluster_LDADD =  ../lib/libfileHandler.la \
		  ../lib/libfileHandler.la \
                  ../lib/libutility.la \
                  ../lib/libmapFactory.la \
                  ../lib/libcsFactory.la \
                  ../lib/libCORRE_MA_OPE.la \
		  ../lib/libimageSource.la\
                  ../lib/libULAGENT.la \
                  ../lib/libsimSystem.la \
		  ../lib/libtimesearch.la\
		  ../lib/libTimeStamp.la \
		  ../lib/libULConstraintSolver.la\
		  ../lib/libPowerUpdater.la\
		  ../lib/libImagePowerUpdater.la\
		  ../lib/libULCS1b.la
