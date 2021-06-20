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
		loadcount = 0;
		dt_start = 0.;
		vtk_timelabel = "TimeValue";
	}
	~VTKreader() {
		delete[] betah;
	}

	
	//vtkSmartPointer<vtkXMLStructuredGridReader> reader0;
	//vtkSmartPointer<vtkXMLStructuredGridReader> reader1;
	vtkXMLStructuredGridReader* reader0;
	vtkXMLStructuredGridReader* reader1;
	//vtkFloatArray* floatdata;
	
	vtkStructuredGrid* dataset0;
	vtkStructuredGrid* dataset1;
	vtkDataArray* U0;
	vtkDataArray* U1;
	vtkCellDataToPointData* c2p0;
	vtkCellDataToPointData* c2p1;

        bool input2d = false;
	int nx, ny, nl;
	double dx, dy;
	double dt_start;
	double bounds[6];
	double* betah; // stretching factor of each cell (percent of total height)
	int switch2d = 0;
	double t0, t1, dt, tmin, tmax, zmin;
	int vfraq_field_located;
	int velo_field_located;
	int dimensions[3];

	int loadcount, Uindex;
	bool cell2Pointdata = false;

	int runOptionSwitch = 0;
	vector<string>* filevec;
	string vtkfilepath, vtk_prefix, Uname, vtk_timelabel;

	double* trilinear_interpolation(double* res, double tpt, double xpt, double ypt, double zpt);

	double* bilinear_interpolation(double* res, double tpt, double xpt, double ypt);

	double* bilinear_interpolation_xy(double* res, double tpt, double xpt, double zpt);
        
    double* linear_interpolation(double* res, double tpt, double xpt);

	bool CheckTime(double tpt);

	void write_vtk(bool endtime);

	void export_vtu(FILE* fp, bool last);

	double z2s(double z, double wave_elev, double depth);
	double s2z(double s, double wave_elev, double depth);

	vector<string>* listdir(const char* dirname, const char* suffix);
	void init(double tpt);
	void update(double tpt);
	double u(double tpt, double xpt, double ypt, double zpt);
	double v(double tpt, double xpt, double ypt, double zpt);
	double w(double tpt, double xpt, double ypt, double zpt);
	double eta(double tpt, double xpt, double ypt);
	double seabed(double xpt, double ypt);
	double getTimeFromVTKFile(string path, const char* fname);
	void loadInit(string path, const char* fname);
	void loadNext(string path, const char* fname);
	double stretchInterpLocatorZ(double x, int* iptr, int nxp, int nyp);
	//void loadNext(string path, const char* fname);


};

/* this class should contain:
1. function which reads the Structured vtk file.
2. function which transforms to sigma space. function should differenciate between vertical lagrangian meshes and regular meshes.
3. locator function (can be taken from lsgrid)
4. interp function (can be taken from lsgrid)

*/
