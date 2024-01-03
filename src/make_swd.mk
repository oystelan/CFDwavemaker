BUILD_DIR += ../builds/linux64/
SWD_INCL = ../swd/inc
CC      := c++
CCFLAGS := -O2 -fPIC -pthread -std=c++11 -fopenmp -DSWD_enable=1 -I$(SWD_INCL)
LDFLAGS := -L./ -L../swd/lib
LIBS += -lm -lgfortran -lfftw3


#export PATH:=$(VTK_INCL):${PATH}

TARGETS:= CFDwavemaker

TARGETS_SHARED_OMP_SWD:= $(addsuffix _swd_openmp.so, $(TARGETS))
TARGETS_STATIC_OMP_SWD:= $(addsuffix _swd_openmp.a, $(TARGETS))

MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ	   := $(MAINS) Stokes5.o Irregular.o Utils.o Wavemaker.o SpectralWaveData.o probes.o lsgrid_spline.o  


.PHONY: clean all

all: $(TARGETS_SHARED_OMP_SWD) $(TARGETS_STATIC_OMP_SWD)

clean:
	rm -f $(OBJ) *f90.o *F90.o

$(OBJ):: %.o : %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS) 

$(TARGETS_SHARED_OMP_SWD): $(OBJ)
	cp ../swd/cpp/SpectralWaveData.cpp .
	cp ../swd/inc/SpectralWaveData.h .
	cp ../swd/inc/spectral_wave_data.h .
	cp ../swd/lib/libSpectralWaveData.a .
	ar x libSpectralWaveData.a
	$(CC) $(CCFLAGS) -shared -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_SWD): $(OBJ)
	rm -f $(BUILD_DIR)lib$@ 
	ar rvs -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o
	chmod 775 $(BUILD_DIR)lib$@
