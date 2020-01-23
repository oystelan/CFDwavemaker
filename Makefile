CC      := g++
CCFLAGS := -O2 -fPIC -std=c++11
LDFLAGS :=
LIBS += -lm
BUILD_DIR += ./builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED:= $(addsuffix .so, $(TARGETS))
TARGETS_STATIC:= $(addsuffix .a, $(TARGETS))
#TARGETS_SHARED_OMP:= $(addsuffix _openmp.so, $(TARGETS))
#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Wavespectra.o Utils.o $(MAINS)
DEPS   := CFDwavemaker.h Stokes5.h Irregular.h Wavespectra.h Utils.h
#OBJ_OMP := $(OBJ)

.PHONY: all clean static shared clean_o

all: $(TARGETS_SHARED) $(TARGETS_STATIC) #clean_o $(TARGETS_SHARED_OMP)

static: $(TARGETS_STATIC)

shared: $(TARGETS_SHARED) #clean_o $(TARGETS_SHARED_OMP)

clean:
	rm -f lib$(TARGETS_SHARED) lib$(TARGETS_STATIC) $(OBJ)

clean_o: rm -f $(OBJ)

$(OBJ): %.o : %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CCFLAGS)

#$(OBJ_OMP): %.o : %.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CCFLAGS) -fopenmp

$(TARGETS_SHARED): $(OBJ)
	$(CC) -shared -o $(BUILD_DIR)lib$@ $(LIBS) $^ $(CCFLAGS_SHARED) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC): $(OBJ)
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@

#$(TARGETS_SHARED_OMP): $(OBJ_OMP)
#	$(CC) -c -o $(OBJ_OMP) $< $(CCFLAGS) -fopenmp
#	$(CC) -shared -fopenmp -o $(BUILD_DIR)lib$@ $(LIBS) $^ $(CCFLAGS_SHARED) $(LDFLAGS)
#	chmod 775 $(BUILD_DIR)lib$@

#$(TARGETS_STATIC_OMP): $(OBJ)
#	ar rvs -o $(BUILD_DIR)lib$@ $^
#	chmod 775 $(BUILD_DIR)lib$@



