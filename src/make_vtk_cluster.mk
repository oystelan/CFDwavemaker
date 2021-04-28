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
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS) -Wl,-rpath $(VTK_LIBS)/libvtkIOImport-9.0.a $(VTK_LIBS)/libvtkFiltersGeneric-9.0.a $(VTK_LIBS)/libvtkIOXML-9.0.a $(VTK_LIBS)/libvtkIOXMLParser-9.0.a $(VTK_LIBS)/libvtkexpat-9.0.a $(VTK_LIBS)/libvtkIOCore-9.0.a $(VTK_LIBS)/libvtklz4-9.0.a $(VTK_LIBS)/libvtklzma-9.0.a $(VTK_LIBS)/libvtklibharu-9.0.a $(VTK_LIBS)/libvtkzlib-9.0.a $(VTK_LIBS)/libvtkFiltersCore-9.0.a $(VTK_LIBS)/libvtkCommonExecutionModel-9.0.a $(VTK_LIBS)/libvtkCommonDataModel-9.0.a $(VTK_LIBS)/libvtkCommonSystem-9.0.a $(VTK_LIBS)/libvtkCommonMisc-9.0.a $(VTK_LIBS)/libvtkCommonTransforms-9.0.a $(VTK_LIBS)/libvtkCommonMath-9.0.a $(VTK_LIBS)/libvtkCommonCore-9.0.a $(VTK_LIBS)/libvtkloguru-9.0.a $(VTK_LIBS)/libvtksys-9.0.a -ldl -lpthread 
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_VTK): EXTRA_FLAGS = -ldl -lpthread -fopenmp
$(TARGETS_STATIC_OMP_VTK): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@
	export PATH=$(PATH):$(VTK_LIBS)
	ar -M <vtk_static_cluster.mri

