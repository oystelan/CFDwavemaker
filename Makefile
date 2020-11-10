CC      := g++
CCFLAGS := -O2 -fPIC -pthread -std=c++11
LDFLAGS :=
LIBS += -lm libSpectralWaveData.a
BUILD_DIR += ./builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED:= $(addsuffix .so, $(TARGETS))
TARGETS_STATIC:= $(addsuffix .a, $(TARGETS))
TARGETS_SHARED_OMP:= $(addsuffix _openmp.so, $(TARGETS))
TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Wavespectra.o Utils.o Wavemaker.o sgrid.o SpectralWaveData.o $(MAINS)
DEPS   := CFDwavemaker.h Stokes5.h Irregular.h Wavespectra.h Utils.h Wavemaker.h sgrid.h SpectralWaveData.h
#OBJ_OMP := $(OBJ)

.PHONY: all clean openmp

all: $(TARGETS_SHARED) $(TARGETS_STATIC)

clean:
	rm -f $(OBJ)

openmp: $(TARGETS_SHARED_OMP) $(TARGETS_STATIC_OMP)


$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

#$(OBJ_OMP): %.o : %.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CCFLAGS) -fopenmp

$(TARGETS_SHARED): $(OBJ)
	$(CC) -shared -o $(BUILD_DIR)lib$@ $(LIBS) $^ $(CCFLAGS_SHARED) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC): $(OBJ)
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_SHARED_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_SHARED_OMP): $(OBJ) 
	$(CC) -shared -fopenmp -o $(BUILD_DIR)lib$@ $(LIBS) $^ $(CCFLAGS_SHARED) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_STATIC_OMP): $(OBJ) 
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@




