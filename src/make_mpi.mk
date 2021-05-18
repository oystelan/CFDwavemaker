CC      := mpic++
CCFLAGS := -O2 -fPIC -pthread -std=c++11
LIBS += -lm 
BUILD_DIR += ../builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED_OMP:= $(addsuffix _mpi.so, $(TARGETS))
TARGETS_STATIC_OMP:= $(addsuffix _mpi.a, $(TARGETS))

MAINS  := $(addsuffix MPI.o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Utils.o Wavemaker.o lsgridMPI.o probes.o $(MAINS) 


.PHONY: clean openmp

all: $(TARGETS_SHARED_OMP) $(TARGETS_STATIC_OMP)

clean:
	rm -f $(OBJ)

openmp: $(TARGETS_SHARED_OMP) $(TARGETS_STATIC_OMP)

$(OBJ):: %.o : %.cpp	 
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

# build with openmp support 
$(TARGETS_SHARED_OMP): $(OBJ)
	$(CC) $(CCFLAGS) -shared -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@
