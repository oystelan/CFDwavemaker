CC      := g++-9
CCFLAGS := -O2 -fPIC -pthread -std=c++11
LDFLAGS := -L./ -L/usr/lib/gcc/x86_64-linux-gnu/9 
LIBS += -lm -lgfortran
BUILD_DIR += ./builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED:= $(addsuffix .so, $(TARGETS))
TARGETS_STATIC:= $(addsuffix .a, $(TARGETS))
TARGETS_SHARED_OMP:= $(addsuffix _openmp.so, $(TARGETS))
TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Wavespectra.o Utils.o Wavemaker.o sgrid.o $(MAINS)
DEPS   := CFDwavemaker.h Stokes5.h Irregular.h Wavespectra.h Utils.h Wavemaker.h sgrid.h
#OBJ_OMP := $(OBJ)

.PHONY: all clean openmp

all: $(TARGETS_SHARED) $(TARGETS_STATIC)

test: $(TARGETS_SHARED)

clean:
	rm -f $(OBJ)
	#ar x libSpectralWaveData.a

openmp: $(TARGETS_SHARED_OMP) $(TARGETS_STATIC_OMP)

static: $(TARGETS_STATIC) $(TARGETS_STATIC_OMP)


$(OBJ): %.o : %.cpp
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

#$(OBJ_OMP): %.o : %.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CCFLAGS) -fopenmp

$(TARGETS_SHARED): $(OBJ) 
	$(CC) $(CCFLAGS) -shared -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC): $(OBJ)
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_SHARED_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_SHARED_OMP): $(OBJ)
	$(CC) $(CCFLAGS) -shared -fopenmp -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP): EXTRA_FLAGS = -fopenmp 
$(TARGETS_STATIC_OMP): $(OBJ) 
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@




