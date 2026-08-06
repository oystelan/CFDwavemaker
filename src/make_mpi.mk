# MPI build (NO OpenMP -- it's either/or: mixing MPI ranks with OpenMP threads
# in the library is asking for trouble).  Distributes the lsGridSpline
# second-order kinematics computation across MPI ranks (each rank computes a
# block, then MPI_Allgatherv).  The host CFD owns MPI_Init/Finalize; this
# library only uses the existing context.  OMPI_SKIP_MPICXX avoids a link
# dependency on OpenMPI's deprecated C++ bindings (-lmpi_cxx).
CC      := mpicxx
CCFLAGS := -O2 -fPIC -pthread -std=c++11 -DMPI_enable=1 -DOMPI_SKIP_MPICXX
LIBS += -lm
BUILD_DIR += ../builds/linux64/

TARGETS:= CFDwavemaker
TARGETS_SHARED_MPI:= $(addsuffix _mpi.so, $(TARGETS))
TARGETS_STATIC_MPI:= $(addsuffix _mpi.a, $(TARGETS))

MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o FentonStream.o Irregular.o Utils.o Wavemaker.o lsgrid_spline.o probes.o $(MAINS)

.PHONY: clean mpi

all: $(TARGETS_SHARED_MPI) $(TARGETS_STATIC_MPI)

clean:
	rm -f $(OBJ)

mpi: $(TARGETS_SHARED_MPI) $(TARGETS_STATIC_MPI)

$(OBJ):: %.o : %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

$(TARGETS_SHARED_MPI): $(OBJ)
	$(CC) $(CCFLAGS) -shared -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_MPI): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@
