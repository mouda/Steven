lib_LTLIBRARIES = libfileHandler.la \
		  libutility.la \
		  libimageSource.la\
		  libmapFactory.la \
		  libimageMapFactory.la\
		  libcsFactory.la \
		  libimageCsFactory.la\
		  libCORRE_MA_OPE.la \
		  libULAGENT.la \
		  libsimSystem.la \
		  libtimesearch.la\
		  libTimeStamp.la \
		  libULConstraintSolver.la\
		  libPowerUpdater.la\
		  libImagePowerUpdater.la\
		  libnewImagePowerUpdater.la\
		  libULCS1b.la
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
	     -lz\
	     -lbz2
libfileHandler_la_SOURCES = fileHandler.h fileHandler.cpp
libfileHandler_la_LIBADD = $(DEPS_LIBS)
libutility_la_SOURCES = utility.h utility.cpp
libutility_la_LIBADD = $(DEPS_LIBS)
libULAGENT_la_SOURCES = ULAGENT.h ULAGENT.cpp
libULAGENT_la_LIBADD = $(DEPS_LIBS)
libsimSystem_la_SOURCES = simSystem.h simSystem.cpp
libsimSystem_la_LIBADD = $(DEPS_LIBS)
libmapFactory_la_SOURCES = map.h map.cpp mapFactory.h mapFactory.cpp utility.h utility.cpp
libmapFactory_la_LIBADD = $(DEPS_LIBS)
libimageMapFactory_la_SOURCES = imageMap.h imageMap.cpp imageMapFactory.h imageMapFactory.cpp utility.h utility.cpp
libimageMapFactory_la_LIBADD = $(DEPS_LIBS)
libcsFactory_la_SOURCES = clusterStructure.h clusterStructure.cpp csFactory.h csFactory.cpp
libcsFactory_la_LIBADD = $(DEPS_LIBS)
libimageCsFactory_la_SOURCES = clusterStructure.h clusterStructure.cpp imageCsFactory.h imageCsFactory.cpp
libimageCsFactory_la_LIBADD = $(DEPS_LIBS) 
libCORRE_MA_OPE_la_SOURCES = CORRE_MA_OPE.h CORRE_MA_OPE.cpp
libCORRE_MA_OPE_la_LIBADD = $(DEPS_LIBS)
libimageSource_la_SOURCES = imageSource.h imageSource.cpp
libimageSource_la_LIBADD = $(DEPS_LIBS) 
libtimesearch_la_SOURCES = timesearch.h timesearch.cpp
libtimesearch_la_LIBADD = $(DEPS_LIBS)
libTimeStamp_la_SOURCES = TimeStamp.h TimeStamp.cpp
libTimeStamp_la_LIBADD = $(DEPS_LIBS)
libULConstraintSolver_la_SOURCES = ULConstraintSolver.h ULConstraintSolver.cpp
libULConstraintSolver_la_LIBADD = $(DEPS_LIBS)
libPowerUpdater_la_SOURCES = powerUpdater.h powerUpdater.cpp
libPowerUpdater_la_LIBADD = $(DEPS_LIBS)
libImagePowerUpdater_la_SOURCES = imagePowerUpdater.h imagePowerUpdater.cpp
libImagePowerUpdater_la_LIBADD = $(DEPS_LIBS)
libnewImagePowerUpdater_la_SOURCES = newImagePowerUpdater.h newImagePowerUpdater.cpp
libnewImagePowerUpdater_la_LIBADD = $(DEPS_LIBS)
libULCS1b_la_SOURCES = ULCS1b.h ULCS1b.cpp
libULCS1b_la_LIBADD = $(DEPS_LIBS)


if DEBUG
AM_CFLAGS = -g3 -O0
AM_CXXFLAGS = -g3 -O0
else
AM_CFLAGS = -O2
AM_CXXFLAGS = -O2
endif
