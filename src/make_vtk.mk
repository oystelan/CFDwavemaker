CC      := c++
CCFLAGS := -O2 -fPIC -pthread -std=c++11 -fopenmp
LDFLAGS := -L./
LIBS += -lm 
BUILD_DIR += ../builds/linux64/
VTK_DIR= /home/oland/vtk/install_vtk-9.0.1
VTK_LIBS = $(VTK_DIR)/lib64
VTK_INCL = $(VTK_DIR)/include/vtk-9.0

export PATH:=$(VTK_INCL):${PATH}

TARGETS:= CFDwavemaker
TARGETS_SHARED_OMP_VTK:= $(addsuffix _vtk_openmp.so, $(TARGETS))
TARGETS_STATIC_OMP_VTK:= $(addsuffix _vtk_openmp.a, $(TARGETS))

#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Utils.o Wavemaker.o lsgrid.o probes.o VTKreader.o lsgrid.o $(MAINS)


.PHONY: all clean 

all: $(TARGETS_SHARED_OMP_VTK) $(TARGETS_STATIC_OMP_VTK)

clean:
	rm -f $(OBJ)

$(OBJ):: %.o : %.cpp	 
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS) -DVTK_enable=1 -I$(VTK_INCL)

$(TARGETS_SHARED_OMP_VTK): $(OBJ)
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS) -ldl -lpthread 
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_VTK): EXTRA_FLAGS = -ldl -lpthread -fopenmp
$(TARGETS_STATIC_OMP_VTK): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@
	ar -M <vtk_static.mri

