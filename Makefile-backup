SUBDIRS = commonLibrary  DirectAcessPowerMax_MC ULSA4b7_DC
export CFLAGS=-O -O1 -O2 -O3

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		make -C $$dir; \
	done; \

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		make  clean -C $$dir; \
	done


