# Compiler
CC      := c++
# VTK input
VTK_VERSION=9.3
VTK_VERSION_LONG=9.3.0
VTK_DIR= /home/oystein/progs/VTK-${VTK_VERSION_LONG}/install
VTK_INCL = $(VTK_DIR)/include/vtk-${VTK_VERSION}
VTK_LIBS = $(VTK_DIR)/lib

SWD_INCL = ../swd/inc

CCFLAGS := -O2 -fPIC -pthread -std=c++11 -fopenmp -DSWD_enable=1 -I$(SWD_INCL) -DVTK_enable=1 -I$(VTK_INCL) 
LDFLAGS := -L./ -L../swd/lib 
LIBS += -lm -lgfortran -lfftw3
BUILD_DIR += ../builds/linux64/
# -D_GLIBCXX_USE_CXX11_ABI=0



export PATH:=$(VTK_INCL):${PATH}

TARGETS:= CFDwavemaker
TARGETS_SHARED_OMP_ALL:= $(addsuffix _all_openmp.so, $(TARGETS))
TARGETS_STATIC_OMP_ALL:= $(addsuffix _all_openmp.a, $(TARGETS))

#TARGETS_STATIC_OMP:= $(addsuffix _openmp.a, $(TARGETS))
MAINS  := $(addsuffix .o, $(TARGETS) )
OBJ    := Stokes5.o Irregular.o Utils.o Wavemaker.o lsgrid.o lsgrid_spline.o SpectralWaveData.o probes.o VTKreader.o $(MAINS)


.PHONY: all clean 

all: $(TARGETS_SHARED_OMP_ALL) $(TARGETS_STATIC_OMP_ALL)

clean:
	rm -f $(OBJ)
	rm -f $(OBJ) $(OBJ_SWD) *f90.o *F90.o
	cp ../swd/cpp/SpectralWaveData.cpp .
	cp ../swd/inc/SpectralWaveData.h .
	cp ../swd/inc/spectral_wave_data.h .
	cp ../swd/lib/libSpectralWaveData.a .

$(OBJ):: %.o : %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) -c -o $@ $< $(CCFLAGS) $(EXTRA_FLAGS)

$(TARGETS_SHARED_OMP_ALL): $(OBJ)
	ar x libSpectralWaveData.a
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o $(LIBS) $(LDFLAGS) -Wl,-rpath $(VTK_LIBS)/libvtkIOImport-${VTK_VERSION}.a $(VTK_LIBS)/libvtkFiltersGeneric-${VTK_VERSION}.a $(VTK_LIBS)/libvtkIOXML-${VTK_VERSION}.a $(VTK_LIBS)/libvtkIOXMLParser-${VTK_VERSION}.a $(VTK_LIBS)/libvtkexpat-${VTK_VERSION}.a $(VTK_LIBS)/libvtkIOCore-${VTK_VERSION}.a $(VTK_LIBS)/libvtklz4-${VTK_VERSION}.a $(VTK_LIBS)/libvtklzma-${VTK_VERSION}.a $(VTK_LIBS)/libvtklibharu-${VTK_VERSION}.a $(VTK_LIBS)/libvtkzlib-${VTK_VERSION}.a $(VTK_LIBS)/libvtkFiltersCore-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonExecutionModel-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonDataModel-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonSystem-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonMisc-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonTransforms-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonMath-${VTK_VERSION}.a $(VTK_LIBS)/libvtkCommonCore-${VTK_VERSION}.a $(VTK_LIBS)/libvtkloguru-${VTK_VERSION}.a $(VTK_LIBS)/libvtksys-${VTK_VERSION}.a $(VTK_LIBS)/libvtkpugixml-${VTK_VERSION}.a -ldl -pthread  
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_ALL): EXTRA_FLAGS = -ldl -pthread -fopenmp
$(TARGETS_STATIC_OMP_ALL): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^ *f90.o *F90.o
	chmod 775 $(BUILD_DIR)lib$@
	export PATH=$(PATH):$(VTK_LIBS)
	echo "open ../builds/linux64/libCFDwavemaker_all_openmp.a" > temp.mri
	echo "addlib ${VTK_LIBS}/libvtkIOCore-${VTK_VERSION}.a" >> temp.mri
	echo "addlib ${VTK_LIBS}/libvtkFiltersCore-${VTK_VERSION}.a" >> temp.mri
	echo "addlib ${VTK_LIBS}/libvtkIOXML-${VTK_VERSION}.a" >> temp.mri
	echo "addlib ${VTK_LIBS}/libvtkIOXMLParser-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtklz4-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtklzma-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkCommonExecutionModel-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkCommonDataModel-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkCommonSystem-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkCommonMisc-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkCommonTransforms-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkCommonMath-${VTK_VERSION}.a" >> temp.mri
	echo "addlib ${VTK_LIBS}/libvtkCommonCore-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkloguru-${VTK_VERSION}.a" >> temp.mri    
	echo "addlib ${VTK_LIBS}/libvtksys-${VTK_VERSION}.a" >> temp.mri  
	echo "addlib ${VTK_LIBS}/libvtkzlib-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkIOImport-${VTK_VERSION}.a" >> temp.mri     
	echo "addlib ${VTK_LIBS}/libvtkFiltersGeneric-${VTK_VERSION}.a" >> temp.mri 
	echo "addlib ${VTK_LIBS}/libvtkFiltersParallel-${VTK_VERSION}.a" >> temp.mri  
	echo "addlib ${VTK_LIBS}/libvtkexpat-${VTK_VERSION}.a" >> temp.mri
	echo "addlib ${VTK_LIBS}/libvtkpugixml-${VTK_VERSION}.a" >> temp.mri
	echo "addlib ${VTK_LIBS}/libvtklibharu-${VTK_VERSION}.a" >> temp.mri    
	echo "save" >> temp.mri
	echo "end" >> temp.mri
	ar -M <temp.mri
	rm temp.mri

