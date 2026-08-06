# Default build: OpenMP + SWD.
# The library version is read from the VERSION file (single source of truth,
# also compiled into the code via -DCFDWM_VERSION). Outputs are versioned
# (libCFDwavemaker_v<X.Y.Z>_openmp.so/.a) with unversioned symlinks
# (libCFDwavemaker_openmp.so/.a) so examples and tests keep working across
# version bumps.
BUILD_DIR += ../builds/linux64/
SWD_INCL = ../swd/inc
VERSION := $(strip $(shell cat VERSION))
CC      := g++
CCFLAGS := -O3 -fPIC -pthread -std=c++11 -fopenmp -DSWD_enable=1 -I$(SWD_INCL) -DCFDWM_VERSION=\"v$(VERSION)\"
LDFLAGS := -L./ -L../swd/lib
LIBS += -lm -lgfortran -lfftw3


#export PATH:=$(VTK_INCL):${PATH}

TARGETS:= CFDwavemaker

TARGETS_SHARED_SWD:= $(addsuffix _v$(VERSION)_openmp.so, $(TARGETS))
TARGETS_STATIC_SWD:= $(addsuffix _v$(VERSION)_openmp.a, $(TARGETS))

MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ	   := $(MAINS) Stokes5.o FentonStream.o Irregular.o Utils.o Wavemaker.o SpectralWaveData.o probes.o lsgrid_spline.o


.PHONY: clean all

all: $(TARGETS_SHARED_SWD) $(TARGETS_STATIC_SWD)

clean:
	rm -f $(OBJ) *f90.o *F90.o

$(OBJ):: %.o : %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS) 

$(TARGETS_SHARED_SWD): $(OBJ)
	cp ../swd/cpp/SpectralWaveData.cpp .
	cp ../swd/inc/SpectralWaveData.h .
	cp ../swd/inc/spectral_wave_data.h .
	cp ../swd/lib/libSpectralWaveData.a .
	ar x libSpectralWaveData.a
	$(CC) $(CCFLAGS) -shared -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@
	ln -sf lib$@ $(BUILD_DIR)libCFDwavemaker_openmp.so

$(TARGETS_STATIC_SWD): $(OBJ)
	rm -f $(BUILD_DIR)lib$@ 
	ar rvs -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o
	chmod 775 $(BUILD_DIR)lib$@
	ln -sf lib$@ $(BUILD_DIR)libCFDwavemaker_openmp.a
