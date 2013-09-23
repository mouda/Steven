SUBDIRS = commonLibrary winULSA2g_MC winULSA2gx_MC winULSA2i2_MC winULSA3g_DC winULSA3gx_DC winULSA3i2_DC winULSA4b2_DC winULSAkmeans2i_MC winULSAkmeans4b2_DC DirectAcessPowerMax_MC
export CFLAGS=-O -O1 -O2 -O3

all: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		make -C $$dir; \
	done; \

clean: $(SUBDIRS)
	for dir in $(SUBDIRS); do \
		make  clean -C $$dir; \
	done


