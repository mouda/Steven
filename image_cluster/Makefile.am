bin_PROGRAMS = genCluster 
INCLUDES =-I../lib
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
		    minPowerImageCluster.h minPowerImageCluster.cpp\
		    minPowerImageCsFactory.h minPowerImageCsFactory.cpp\
		    baselineImageCluster.h baselineImageCluster.cpp\
		    baselineImageCsFactory.h baselineImageCsFactory.cpp\
		    simulator.h simulator.cpp 

LOCAL_INCLUDE = -I../lib/
if DEBUG
AM_CFLAGS = $(IPOPT_CPPFLAGS) $(IPOPT_LIBS) $(LOCAL_INCLUDE)  -I../ -g3 -ggdb -O0
AM_CXXFLAGS = $(IPOPT_CPPFLAGS) $(IPOPT_LIBS) $(LOCAL_INCLUDE) -I../ -g3 -ggdb -O0
else
AM_CFLAGS = $(IPOPT_CPPFLAGS) $(IPOPT_LIBS) $(LOCAL_INCLUDE)  -I../ -O2
AM_CXXFLAGS = $(IPOPT_CPPFLAGS) $(IPOPT_LIBS) $(LOCAL_INCLUDE)  -I../ -O2
endif
genCluster_LDADD = $(IPOPT_CPPFLAGS) $(IPOPT_LIBS) ../lib/libfileHandler.la \
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
		  ../lib/libimageMapFactory.la\
		  ../lib/libnewImagePowerUpdater.la\
		  ../lib/libimageCsFactory.la\
		  ../lib/libULCS1b.la
