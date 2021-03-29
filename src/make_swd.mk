BUILD_DIR += ../builds/linux64/
SWD_INCL = ../swd/inc
CC      := g++-9
CCFLAGS := -O2 -fPIC -pthread -std=c++17 -fopenmp -DSWD_enable=1 -I$(SWD_INCL)
LDFLAGS := -L./ -L/home/oland/programs/CFDwavemaker/swd/lib -L/usr/lib/gcc/x86_64-linux-gnu/9
LIBS += -lm -lgfortran


export PATH:=$(VTK_INCL):${PATH}

TARGETS:= CFDwavemaker

TARGETS_SHARED_OMP_SWD:= $(addsuffix _swd_openmp.so, $(TARGETS))
TARGETS_STATIC_OMP_SWD:= $(addsuffix _swd_openmp.a, $(TARGETS))

MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ	   := $(MAINS) Stokes5.o Irregular.o Utils.o Wavemaker.o SpectralWaveData.o probes.o lsgrid.o  


.PHONY: clean all

all: $(TARGETS_SHARED_OMP_SWD) $(TARGETS_STATIC_OMP_SWD)

clean:
	rm -f $(OBJ) $(OBJ_SWD) *f90.o *F90.o
	cp ../swd/cpp/SpectralWaveData.cpp .
	cp ../swd/inc/SpectralWaveData.h .
	cp ../swd/inc/spectral_wave_data.h .

$(OBJ):: %.o : %.cpp	 
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS) 

$(TARGETS_SHARED_OMP_SWD): $(OBJ) 
	ar x libSpectralWaveData.a
	$(CC) $(CCFLAGS) -shared -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_SWD): $(OBJ) 
	ar rvs -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o
	chmod 775 $(BUILD_DIR)lib$@
