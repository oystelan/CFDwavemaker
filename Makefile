CC      := g++-9
CCFLAGS := -O2 -fPIC -pthread -std=c++17 -DSWD_enable=1
LDFLAGS := -L./ -L./swd/lib -L/usr/lib/gcc/x86_64-linux-gnu/9 
LIBS += -lm -lgfortran
BUILD_DIR += ./builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED:= $(addsuffix .so, $(TARGETS))
TARGETS_STATIC:= $(addsuffix .a, $(TARGETS))
TARGETS_SHARED_OMP:= $(addsuffix _openmp.so, $(TARGETS))
TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Utils.o Wavemaker.o lsgrid.o SpectralWaveData.o probes.o $(MAINS)
DEPS   := CFDwavemaker.h Stokes5.h Irregular.h Utils.h Wavemaker.h lsgrid.h probes.h spectral_wave_data.h SpectralWaveData.h
#OBJ_OMP := $(OBJ)

.PHONY: all clean openmp

all: $(TARGETS_SHARED) $(TARGETS_STATIC)

test: $(TARGETS_SHARED)

clean:
	rm -f $(OBJ)
	ar x ./swd/lib/libSpectralWaveData.a
	cp ./swd/cpp/SpectralWaveData.cpp .
	cp ./swd/inc/SpectralWaveData.h .
	cp ./swd/inc/spectral_wave_data.h .

openmp: $(TARGETS_SHARED_OMP) $(TARGETS_STATIC_OMP)

static: $(TARGETS_STATIC) $(TARGETS_STATIC_OMP)


$(OBJ): %.o : %.cpp
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

#$(OBJ_OMP): %.o : %.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CCFLAGS) -fopenmp

$(TARGETS_SHARED): $(OBJ) *.o
	$(CC) $(CCFLAGS) -shared -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC): $(OBJ) *.o
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_SHARED_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_SHARED_OMP): $(OBJ) *.o
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_STATIC_OMP): $(OBJ) *.o
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@




