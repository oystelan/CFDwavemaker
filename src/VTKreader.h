#pragma once

// VTK includes
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellDataToPointData.h>

// Other includes for system access
#include <vector>
#include <string>
#include <cstring>
//#include <algorithm>
//#include <sys/types.h>
#include "dirent.h"
#include <iterator>
//#include <cstddef>

using namespace std;

class VTKreader
{
private:
	


public:
	VTKreader() {
		Uindex = -1;
	}
	~VTKreader() {
		delete[] ZB;
	}

	
	//vtkSmartPointer<vtkXMLStructuredGridReader> reader0;
	//vtkSmartPointer<vtkXMLStructuredGridReader> reader1;
	vtkXMLStructuredGridReader* reader0;
	vtkXMLStructuredGridReader* reader1;
	vtkFloatArray* floatdata;
	vtkDoubleArray* doubledata;
	vtkStructuredGrid* dataset0;
	vtkStructuredGrid* dataset1;

	int nx, ny, nl, dx, dy;
	double bounds[6];
	double* ZB;
	double* beta; // stretching factor of each cell (percent of total height)
	int switch2d = 0;
	double t0, t1, dt;
	int vfraq_field_located;
	int velo_field_located;
	int dimensions[3];

	int filecount, Uindex;
	bool cell2Pointdata = false;

	int runOptionSwitch = 0;
	vector<string>* filevec;
	string vtkfilepath, vtk_prefix, Uname;

	double* trilinear_interpolation(double tpt, double xpt, double ypt, double zpt);



	double z2s(double z, double wave_elev, double depth);

	double s2z(double s, double wave_elev, double depth);


	vector<string>* listdir(const char* dirname, const char* suffix);
	void init();
	void loadInit(string path, const char* fname);
	//void loadNext(string path, const char* fname);


};

/* this class should contain:
1. function which reads the Structured vtk file.
2. function which transforms to sigma space. function should differenciate between vertical lagrangian meshes and regular meshes.
3. locator function (can be taken from lsgrid)
4. interp function (can be taken from lsgrid)

*/