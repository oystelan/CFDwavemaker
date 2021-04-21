CC      := c++
CCFLAGS := -O2 -fPIC -pthread -std=c++11 -fopenmp
LDFLAGS := -L./ -L/usr/lib/gcc/x86_64-linux-gnu/9
LIBS += -lm 
BUILD_DIR += ../builds/linux64/
VTK_DIR= /home/oland/programs/vtk-9.0.1/installdir
VTK_LIBS = $(VTK_DIR)/lib
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
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS) -Wl,-rpath,/usr/lib/nvidia-384: /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkWrappingTools-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkViewsInfovis-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonColor-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkViewsContext2D-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkloguru-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkTestingRendering-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtksys-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingVolumeOpenGL2-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingOpenGL2-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkglew-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingLabel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingLOD-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingImage-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOVeraOut-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOTecplotTable-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOSegY-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOParallelXML-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOPLY-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOOggTheora-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtktheora-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkogg-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIONetCDF-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtknetcdf-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOMotionFX-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOParallel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkjsoncpp-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOMINC-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOLSDyna-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOInfovis-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtklibxml2-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkzlib-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOImport-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOGeometry-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOVideo-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOMovie-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOExportPDF-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtklibharu-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOExportGL2PS-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingGL2PSOpenGL2-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkgl2ps-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkpng-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOExport-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingVtkJS-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingSceneGraph-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOExodus-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkexodusII-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOEnSight-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOCityGML-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkpugixml-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOAsynchronous-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOAMR-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkInteractionImage-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingStencil-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingStatistics-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingMorphological-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingMath-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOSQL-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtksqlite-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkGeovisCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtklibproj-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkInfovisLayout-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkViewsCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkInteractionWidgets-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingVolume-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingAnnotation-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingHybrid-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingColor-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkInteractionStyle-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersTopology-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersSelection-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersSMP-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersProgrammable-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersPoints-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersVerdict-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkverdict-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersParallelImaging-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersImaging-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingGeneral-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersHyperTree-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersGeneric-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersFlowPaths-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersAMR-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersParallel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersTexture-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersModeling-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersHybrid-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingUI-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkDomainsChemistry-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkChartsCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkInfovisCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersExtraction-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkParallelDIY-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOXML-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOXMLParser-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkexpat-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkParallelCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOLegacy-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkdoubleconversion-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtklz4-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtklzma-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersStatistics-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingFourier-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingSources-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkIOImage-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkDICOMParser-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkjpeg-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkmetaio-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtktiff-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingContext2D-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingFreeType-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkfreetype-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersSources-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkImagingCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersGeometry-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersGeneral-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonComputationalGeometry-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonExecutionModel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonDataModel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonSystem-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonMisc-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonTransforms-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonMath-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtklibharu-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingOpenGL2-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkglew-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkjsoncpp-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtknetcdf-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkhdf5_hl-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkhdf5-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingUI-9.0.a /usr/lib/x86_64-linux-gnu/libX11.so /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkpng-9.0.a -lm /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkpugixml-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkjpeg-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkzlib-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkRenderingCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonColor-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersGeometry-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersSources-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersGeneral-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonComputationalGeometry-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkFiltersCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonExecutionModel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonDataModel-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonSystem-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonMisc-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonTransforms-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonMath-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkCommonCore-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtkloguru-9.0.a /home/oland/programs/vtk-9.0.1/build_VTK-9.0.1/lib/libvtksys-9.0.a -ldl -lpthread 
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_VTK): EXTRA_FLAGS = -ldl -lpthread -fopenmp
$(TARGETS_STATIC_OMP_VTK): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@
	ar -M <vtk_static.mri

