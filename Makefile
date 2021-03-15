CC      := g++-9
CCFLAGS := -O2 -fPIC -pthread -std=c++17
LDFLAGS := -L./ -L./swd/lib -L/usr/lib/gcc/x86_64-linux-gnu/9 
LIBS += -lm -lgfortran
BUILD_DIR += ./builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED:= $(addsuffix .so, $(TARGETS))
TARGETS_STATIC:= $(addsuffix .a, $(TARGETS))
TARGETS_SHARED_OMP:= $(addsuffix _openmp.so, $(TARGETS))
TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
TARGETS_SHARED_OMP_SWD:= $(addsuffix _swd_openmp.so, $(TARGETS))
TARGETS_STATIC_OMP_SWD:= $(addsuffix _swd_openmp.a, $(TARGETS))

#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Utils.o Wavemaker.o lsgrid.o probes.o $(MAINS)
OBJ_HOSM := SpectralWaveData.o
DEPS   := CFDwavemaker.h Stokes5.h Irregular.h Utils.h Wavemaker.h lsgrid.h probes.h 
DEPS_HOSM := spectral_wave_data.h SpectralWaveData.h

.PHONY: all clean openmp openmp_swd

all: $(TARGETS_SHARED) $(TARGETS_STATIC)

test: $(TARGETS_SHARED)

clean:
	rm -f $(OBJ)
	ar x ./swd/lib/libSpectralWaveData.a
	cp ./swd/cpp/SpectralWaveData.cpp .
	cp ./swd/inc/SpectralWaveData.h .
	cp ./swd/inc/spectral_wave_data.h .

openmp: $(TARGETS_SHARED_OMP) $(TARGETS_STATIC_OMP)

openmp_swd: $(TARGETS_SHARED_OMP_SWD) $(TARGETS_STATIC_OMP_SWD)

static: $(TARGETS_STATIC) $(TARGETS_STATIC_OMP)


$(OBJ): %.o : %.cpp
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

$(OBJ_HOSM): %.o : %.cpp
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

# build serial
$(TARGETS_SHARED): $(OBJ) *.o
	$(CC) $(CCFLAGS) -shared -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC): $(OBJ) *.o
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@


# build with openmp support
$(TARGETS_SHARED_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_SHARED_OMP): $(OBJ) 
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_STATIC_OMP): $(OBJ) 
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@

# Build CFDwavemaker with openmp and including SWD library
$(TARGETS_SHARED_OMP_SWD): EXTRA_FLAGS = -fopenmp -DSWD_enable=1
$(TARGETS_SHARED_OMP_SWD): $(OBJ) $(OBJ_SWD) *f90.o
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_SWD): EXTRA_FLAGS = -fopenmp -DSWD_enable=1
$(TARGETS_STATIC_OMP_SWD): $(OBJ) $(OBJ_SWD) *f90.o
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@

