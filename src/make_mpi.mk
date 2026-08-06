# MPI build (NO OpenMP -- it's either/or: mixing MPI ranks with OpenMP threads
# in the library is asking for trouble).  Distributes the lsGridSpline
# second-order kinematics computation across MPI ranks (each rank computes a
# block, then MPI_Allgatherv).  The host CFD owns MPI_Init/Finalize; this
# library only uses the existing context.  OMPI_SKIP_MPICXX avoids a link
# dependency on OpenMPI's deprecated C++ bindings (-lmpi_cxx).
#
# SWD is included (like the default build): each MPI rank is a separate
# process with its own swd reader object, so the swd+lsgrid path (wavetype
# 34) fills its column block per rank and assembles with MPI_Allgatherv.
BUILD_DIR += ../builds/linux64/
SWD_INCL = ../swd/inc
CC      := mpicxx
CCFLAGS := -O2 -fPIC -pthread -std=c++11 -DMPI_enable=1 -DOMPI_SKIP_MPICXX -DSWD_enable=1 -I$(SWD_INCL)
LDFLAGS := -L./ -L../swd/lib
LIBS += -lm -lgfortran -lfftw3

TARGETS:= CFDwavemaker
TARGETS_SHARED_MPI:= $(addsuffix _mpi.so, $(TARGETS))
TARGETS_STATIC_MPI:= $(addsuffix _mpi.a, $(TARGETS))

MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o FentonStream.o Irregular.o Utils.o Wavemaker.o SpectralWaveData.o lsgrid_spline.o probes.o $(MAINS)

.PHONY: clean mpi

all: $(TARGETS_SHARED_MPI) $(TARGETS_STATIC_MPI)

clean:
	rm -f $(OBJ) *f90.o *F90.o

mpi: $(TARGETS_SHARED_MPI) $(TARGETS_STATIC_MPI)

$(OBJ):: %.o : %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

$(TARGETS_SHARED_MPI): $(OBJ)
	cp ../swd/cpp/SpectralWaveData.cpp .
	cp ../swd/inc/SpectralWaveData.h .
	cp ../swd/inc/spectral_wave_data.h .
	cp ../swd/lib/libSpectralWaveData.a .
	ar x libSpectralWaveData.a
	$(CC) $(CCFLAGS) -shared -fPIC -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o $(LIBS) $(LDFLAGS)
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_MPI): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o
	chmod 775 $(BUILD_DIR)lib$@
