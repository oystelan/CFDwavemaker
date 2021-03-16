#pragma once

// VTK includes
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>

// Other includes for system access
#include <vector>
#include <string>
#include <cstring>
//#include <algorithm>
//#include <sys/types.h>
#include "dirent.h"
//#include <iterator>
//#include <cstddef>

using namespace std;

class VTKreader
{
private:


public:
	vtkSmartPointer<vtkXMLStructuredGridReader> reader0;
	vtkSmartPointer<vtkXMLStructuredGridReader> reader1;
	vtkFloatArray* floatdata;
	vtkDoubleArray* doubledata;
	vtkStructuredGrid* dataset0;
	vtkStructuredGrid* dataset1;

	double bounds[6];
	int switch2d = 0;
	double time0, time1;
	int vfraq_field_located;
	int velo_field_located;

	int runOptionSwitch = 0;
	vector<string>* filevec;
	string vtkfilepath, vtk_prefix, chF, chV, chP;

	vector<string>* listdir(const char* dirname, const char* suffix, int& IERR);
	void loadInit(string path, const char* fname);
	void loadNext(string path, const char* fname);


};

/* this class should contain:
1. function which reads the Structured vtk file.
2. function which transforms to sigma space. function should differenciate between vertical lagrangian meshes and regular meshes.
3. locator function (can be taken from lsgrid)
4. interp function (can be taken from lsgrid)

*/