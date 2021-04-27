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
	$(CC) $(CCFLAGS) -shared -fopenmp -fPIC -o $(BUILD_DIR)lib$@ $^ $(LIBS) $(LDFLAGS) -Wl,-rpath $(VTK_LIBS)/libvtkWrappingTools-9.0.a $(VTK_LIBS)/libvtkViewsInfovis-9.0.a $(VTK_LIBS)/libvtkCommonColor-9.0.a $(VTK_LIBS)/libvtkViewsContext2D-9.0.a $(VTK_LIBS)/libvtkloguru-9.0.a $(VTK_LIBS)/libvtkTestingRendering-9.0.a $(VTK_LIBS)/libvtksys-9.0.a $(VTK_LIBS)/libvtkRenderingVolumeOpenGL2-9.0.a $(VTK_LIBS)/libvtkRenderingOpenGL2-9.0.a $(VTK_LIBS)/libvtkglew-9.0.a $(VTK_LIBS)/libvtkRenderingLabel-9.0.a $(VTK_LIBS)/libvtkRenderingLOD-9.0.a $(VTK_LIBS)/libvtkRenderingImage-9.0.a $(VTK_LIBS)/libvtkIOVeraOut-9.0.a $(VTK_LIBS)/libvtkIOTecplotTable-9.0.a $(VTK_LIBS)/libvtkIOSegY-9.0.a $(VTK_LIBS)/libvtkIOParallelXML-9.0.a $(VTK_LIBS)/libvtkIOPLY-9.0.a $(VTK_LIBS)/libvtkIOOggTheora-9.0.a $(VTK_LIBS)/libvtktheora-9.0.a $(VTK_LIBS)/libvtkogg-9.0.a $(VTK_LIBS)/libvtkIONetCDF-9.0.a $(VTK_LIBS)/libvtknetcdf-9.0.a $(VTK_LIBS)/libvtkIOMotionFX-9.0.a $(VTK_LIBS)/libvtkIOParallel-9.0.a $(VTK_LIBS)/libvtkjsoncpp-9.0.a $(VTK_LIBS)/libvtkIOMINC-9.0.a $(VTK_LIBS)/libvtkIOLSDyna-9.0.a $(VTK_LIBS)/libvtkIOInfovis-9.0.a $(VTK_LIBS)/libvtklibxml2-9.0.a $(VTK_LIBS)/libvtkzlib-9.0.a $(VTK_LIBS)/libvtkIOImport-9.0.a $(VTK_LIBS)/libvtkIOGeometry-9.0.a $(VTK_LIBS)/libvtkIOVideo-9.0.a $(VTK_LIBS)/libvtkIOMovie-9.0.a $(VTK_LIBS)/libvtkIOExportPDF-9.0.a $(VTK_LIBS)/libvtklibharu-9.0.a $(VTK_LIBS)/libvtkIOExportGL2PS-9.0.a $(VTK_LIBS)/libvtkRenderingGL2PSOpenGL2-9.0.a $(VTK_LIBS)/libvtkgl2ps-9.0.a $(VTK_LIBS)/libvtkpng-9.0.a $(VTK_LIBS)/libvtkIOExport-9.0.a $(VTK_LIBS)/libvtkRenderingVtkJS-9.0.a $(VTK_LIBS)/libvtkRenderingSceneGraph-9.0.a $(VTK_LIBS)/libvtkIOExodus-9.0.a $(VTK_LIBS)/libvtkexodusII-9.0.a $(VTK_LIBS)/libvtkIOEnSight-9.0.a $(VTK_LIBS)/libvtkIOCityGML-9.0.a $(VTK_LIBS)/libvtkpugixml-9.0.a $(VTK_LIBS)/libvtkIOAsynchronous-9.0.a $(VTK_LIBS)/libvtkIOAMR-9.0.a $(VTK_LIBS)/libvtkInteractionImage-9.0.a $(VTK_LIBS)/libvtkImagingStencil-9.0.a $(VTK_LIBS)/libvtkImagingStatistics-9.0.a $(VTK_LIBS)/libvtkImagingMorphological-9.0.a $(VTK_LIBS)/libvtkImagingMath-9.0.a $(VTK_LIBS)/libvtkIOSQL-9.0.a $(VTK_LIBS)/libvtksqlite-9.0.a $(VTK_LIBS)/libvtkGeovisCore-9.0.a $(VTK_LIBS)/libvtklibproj-9.0.a $(VTK_LIBS)/libvtkInfovisLayout-9.0.a $(VTK_LIBS)/libvtkViewsCore-9.0.a $(VTK_LIBS)/libvtkInteractionWidgets-9.0.a $(VTK_LIBS)/libvtkRenderingVolume-9.0.a $(VTK_LIBS)/libvtkRenderingAnnotation-9.0.a $(VTK_LIBS)/libvtkImagingHybrid-9.0.a $(VTK_LIBS)/libvtkImagingColor-9.0.a $(VTK_LIBS)/libvtkInteractionStyle-9.0.a $(VTK_LIBS)/libvtkFiltersTopology-9.0.a $(VTK_LIBS)/libvtkFiltersSelection-9.0.a $(VTK_LIBS)/libvtkFiltersSMP-9.0.a $(VTK_LIBS)/libvtkFiltersProgrammable-9.0.a $(VTK_LIBS)/libvtkFiltersPoints-9.0.a $(VTK_LIBS)/libvtkFiltersVerdict-9.0.a $(VTK_LIBS)/libvtkverdict-9.0.a $(VTK_LIBS)/libvtkFiltersParallelImaging-9.0.a $(VTK_LIBS)/libvtkFiltersImaging-9.0.a $(VTK_LIBS)/libvtkImagingGeneral-9.0.a $(VTK_LIBS)/libvtkFiltersHyperTree-9.0.a $(VTK_LIBS)/libvtkFiltersGeneric-9.0.a $(VTK_LIBS)/libvtkFiltersFlowPaths-9.0.a $(VTK_LIBS)/libvtkFiltersAMR-9.0.a $(VTK_LIBS)/libvtkFiltersParallel-9.0.a $(VTK_LIBS)/libvtkFiltersTexture-9.0.a $(VTK_LIBS)/libvtkFiltersModeling-9.0.a $(VTK_LIBS)/libvtkFiltersHybrid-9.0.a $(VTK_LIBS)/libvtkRenderingUI-9.0.a $(VTK_LIBS)/libvtkDomainsChemistry-9.0.a $(VTK_LIBS)/libvtkChartsCore-9.0.a $(VTK_LIBS)/libvtkInfovisCore-9.0.a $(VTK_LIBS)/libvtkFiltersExtraction-9.0.a $(VTK_LIBS)/libvtkParallelDIY-9.0.a $(VTK_LIBS)/libvtkIOXML-9.0.a $(VTK_LIBS)/libvtkIOXMLParser-9.0.a $(VTK_LIBS)/libvtkexpat-9.0.a $(VTK_LIBS)/libvtkParallelCore-9.0.a $(VTK_LIBS)/libvtkIOLegacy-9.0.a $(VTK_LIBS)/libvtkIOCore-9.0.a $(VTK_LIBS)/libvtkdoubleconversion-9.0.a $(VTK_LIBS)/libvtklz4-9.0.a $(VTK_LIBS)/libvtklzma-9.0.a $(VTK_LIBS)/libvtkFiltersStatistics-9.0.a $(VTK_LIBS)/libvtkImagingFourier-9.0.a $(VTK_LIBS)/libvtkImagingSources-9.0.a $(VTK_LIBS)/libvtkIOImage-9.0.a $(VTK_LIBS)/libvtkDICOMParser-9.0.a $(VTK_LIBS)/libvtkjpeg-9.0.a $(VTK_LIBS)/libvtkmetaio-9.0.a $(VTK_LIBS)/libvtktiff-9.0.a $(VTK_LIBS)/libvtkRenderingContext2D-9.0.a $(VTK_LIBS)/libvtkRenderingFreeType-9.0.a $(VTK_LIBS)/libvtkfreetype-9.0.a $(VTK_LIBS)/libvtkRenderingCore-9.0.a $(VTK_LIBS)/libvtkFiltersSources-9.0.a $(VTK_LIBS)/libvtkImagingCore-9.0.a $(VTK_LIBS)/libvtkFiltersGeometry-9.0.a $(VTK_LIBS)/libvtkFiltersGeneral-9.0.a $(VTK_LIBS)/libvtkCommonComputationalGeometry-9.0.a $(VTK_LIBS)/libvtkFiltersCore-9.0.a $(VTK_LIBS)/libvtkCommonExecutionModel-9.0.a $(VTK_LIBS)/libvtkCommonDataModel-9.0.a $(VTK_LIBS)/libvtkCommonSystem-9.0.a $(VTK_LIBS)/libvtkCommonMisc-9.0.a $(VTK_LIBS)/libvtkCommonTransforms-9.0.a $(VTK_LIBS)/libvtkCommonMath-9.0.a $(VTK_LIBS)/libvtkCommonCore-9.0.a $(VTK_LIBS)/libvtklibharu-9.0.a $(VTK_LIBS)/libvtkRenderingOpenGL2-9.0.a $(VTK_LIBS)/libvtkglew-9.0.a $(VTK_LIBS)/libvtkjsoncpp-9.0.a $(VTK_LIBS)/libvtknetcdf-9.0.a $(VTK_LIBS)/libvtkhdf5_hl-9.0.a $(VTK_LIBS)/libvtkhdf5-9.0.a $(VTK_LIBS)/libvtkRenderingUI-9.0.a $(VTK_LIBS)/libvtkpng-9.0.a -lm $(VTK_LIBS)/libvtkpugixml-9.0.a $(VTK_LIBS)/libvtkjpeg-9.0.a $(VTK_LIBS)/libvtkzlib-9.0.a $(VTK_LIBS)/libvtkRenderingCore-9.0.a $(VTK_LIBS)/libvtkCommonColor-9.0.a $(VTK_LIBS)/libvtkFiltersGeometry-9.0.a $(VTK_LIBS)/libvtkFiltersSources-9.0.a $(VTK_LIBS)/libvtkFiltersGeneral-9.0.a $(VTK_LIBS)/libvtkCommonComputationalGeometry-9.0.a $(VTK_LIBS)/libvtkFiltersCore-9.0.a $(VTK_LIBS)/libvtkCommonExecutionModel-9.0.a $(VTK_LIBS)/libvtkCommonDataModel-9.0.a $(VTK_LIBS)/libvtkCommonSystem-9.0.a $(VTK_LIBS)/libvtkCommonMisc-9.0.a $(VTK_LIBS)/libvtkCommonTransforms-9.0.a $(VTK_LIBS)/libvtkCommonMath-9.0.a $(VTK_LIBS)/libvtkCommonCore-9.0.a $(VTK_LIBS)/libvtkloguru-9.0.a $(VTK_LIBS)/libvtksys-9.0.a -ldl -lpthread 
	chmod 775 $(BUILD_DIR)lib$@

$(TARGETS_STATIC_OMP_VTK): EXTRA_FLAGS = -ldl -lpthread -fopenmp
$(TARGETS_STATIC_OMP_VTK): $(OBJ)
	rm -f $(BUILD_DIR)lib$@
	ar rvs -o $(BUILD_DIR)lib$@ $^
	chmod 775 $(BUILD_DIR)lib$@
	export PATH=$(PATH):$(VTK_LIBS)
	ar -M <vtk_static_cluster.mri

