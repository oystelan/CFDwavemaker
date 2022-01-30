CC      := c++
VTK_DIR= /home/oland/programs/vtk/install_vtk9.1.0
SWD_INCL = ../swd/inc
VTK_INCL = $(VTK_DIR)/include/vtk-9.1
CCFLAGS := -O2 -fPIC -pthread -std=c++11 -fopenmp -DSWD_enable=1 -I$(SWD_INCL) -DVTK_enable=1 -I$(VTK_INCL) 
LDFLAGS := -L./ -L../swd/lib 
LIBS += -lm -lgfortran -lfftw3
BUILD_DIR += ../builds/linux64/
# -D_GLIBCXX_USE_CXX11_ABI=0
VTK_LIBS = $(VTK_DIR)/lib64


export PATH:=$(VTK_INCL):${PATH}

TARGETS:= CFDwavemaker
TARGETS_SHARED_OMP_ALL:= $(addsuffix _all_openmp.so, $(TARGETS))
TARGETS_STATIC_OMP_ALL:= $(addsuffix _all_openmp.a, $(TARGETS))

#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Utils.o Wavemaker.o lsgrid.o SpectralWaveData.o probes.o VTKreader.o $(MAINS)


.PHONY: all clean 

all: $(TARGETS_SHARED_OMP_ALL) $(TARGETS_STATIC_OMP_ALL)

clean:
	rm -f $(OBJ)
	rm -f $(OBJ) $(OBJ_SWD) *f90.o *F90.o
	cp ../swd/cpp/SpectralWaveData.cpp .
	cp ../swd/inc/SpectralWaveData.h .
	cp ../swd/inc/spectral_wave_data.h .

$(OBJ):: %.o : %.cpp	 
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

$(TARGETS_SHARED_OMP_ALL): $(OBJ)
	ar x libSpectralWaveData.a
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o $(LIBS) $(LDFLAGS) -Wl,-rpath $(VTK_LIBS)/libvtkIOImport-9.1.a $(VTK_LIBS)/libvtkFiltersGeneric-9.1.a $(VTK_LIBS)/libvtkIOXML-9.1.a $(VTK_LIBS)/libvtkIOXMLParser-9.1.a $(VTK_LIBS)/libvtkexpat-9.1.a $(VTK_LIBS)/libvtkIOCore-9.1.a $(VTK_LIBS)/libvtklz4-9.1.a $(VTK_LIBS)/libvtklzma-9.1.a $(VTK_LIBS)/libvtklibharu-9.1.a $(VTK_LIBS)/libvtkzlib-9.1.a $(VTK_LIBS)/libvtkFiltersCore-9.1.a $(VTK_LIBS)/libvtkCommonExecutionModel-9.1.a $(VTK_LIBS)/libvtkCommonDataModel-9.1.a $(VTK_LIBS)/libvtkCommonSystem-9.1.a $(VTK_LIBS)/libvtkCommonMisc-9.1.a $(VTK_LIBS)/libvtkCommonTransforms-9.1.a $(VTK_LIBS)/libvtkCommonMath-9.1.a $(VTK_LIBS)/libvtkCommonCore-9.1.a $(VTK_LIBS)/libvtkloguru-9.1.a $(VTK_LIBS)/libvtksys-9.1.a $(VTK_LIBS)/libvtkpugixml-9.1.a -ldl -pthread  
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_ALL): EXTRA_FLAGS = -ldl -pthread -fopenmp
$(TARGETS_STATIC_OMP_ALL): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o
	chmod 775 $(BUILD_DIR)lib$@
	export PATH=$(PATH):$(VTK_LIBS)
	ar -M <all_static_cluster.mri

