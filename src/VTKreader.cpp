#include "VTKreader.h"
#include "Utils.h"
#include <sys/types.h>
#include <sys/stat.h>

// Transforms normal z axis to streched sigma coordinates 
// s defined between -1 (seabed) and 0 (free surface)
double VTKreader::z2s(double z, double wave_elev, double depth) {
	return (z - wave_elev) / (depth + wave_elev);
}

// Transforms stretched sigma coordinate to normal z
double VTKreader::s2z(double s, double wave_elev, double depth) {
	return wave_elev + s * (depth + wave_elev);
}


// Function for listing vtk files in directory
vector<string>* VTKreader::listdir(const char* dirname, const char* suffix) {
	DIR* dp;
	dirent* d;
	vector<string>* vec = new vector<string>;

	dp = opendir(dirname);
	if (dp) {
		while ((d = readdir(dp)) != NULL) {

			// you want here. Something like:
			if (strstr(d->d_name, suffix)) {
				//cout << "found an .cxx file: " << d->d_name << "\n";
				vec->push_back(d->d_name);
			}
		}
		sort(vec->begin(), vec->end());
	}
	else {
		cerr << "could not find any files in specified path..." << endl;
		exit(-1);
	}
	return vec;
}

// VTK reader
void VTKreader::init() {

	filevec = listdir(vtkfilepath.c_str(), vtk_prefix.c_str());
	ostream_iterator<string> out(cout, " ");
	cout << "The following list of files identified in specified folder:" << endl;
	copy(filevec->begin(), filevec->end(), out);
	cout << endl;

	//cout << filevec->size();
	if (filevec->size() >= 2) {
		// load vtufiles
		loadInit(vtkfilepath, filevec->at(0).c_str());
		loadNext(vtkfilepath, filevec->at(1).c_str());
		write_vtk(false);
		cout << "VTU Files found and loaded..." << endl;
	}
	else {
		cout << "Not enough (or no) VTU files in folder. A minimum of two files are required.\n\
			Perhaps the wrong path is specified? just trying to be helpful..." << endl;
	}

}

// Function loading data from the frist file into data set 2. This function is only called during initialization
void VTKreader::loadInit(string path, const char* fname) {
	// Load file 1

	path.append(fname);
	cout << path << endl;

	reader1 = vtkXMLStructuredGridReader::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();
	dataset1 = reader1->GetOutput();

	// Get time
	vtkDoubleArray* doubledata = vtkDoubleArray::SafeDownCast(dataset1->GetFieldData()->GetArray("TimeValue"));
	t1 = doubledata->GetValue(0);
	dt = t1 - t0;


	// get bounds
	dataset1->GetBounds(bounds);
	cout << "The bounds: " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << " " << endl;
	if (bounds[2] == bounds[3]) {
		cout << "2D modus switched on." << endl;
		switch2d = 1;
	}

	//get extent of grid (ncells in each direction
	int* dims = dataset1->GetDimensions();
	cout << "Dimensions of grid: " << dims[0] << ", " << dims[1] << ", " << dims[2] << ", " << endl;
	nx = dims[0];
	ny = dims[1];
	nl = dims[2];

	// find cell size dx and dy
	double coord[3];
	dataset1->GetPoint(0, 0, 0, coord);
	double x0 = coord[0];
	double y0 = coord[1];
	cout << coord[0] << " " << coord[1] << " " << coord[2] << endl;
	dataset1->GetPoint(1, 1, 0, coord);
	cout << coord[0] << " " << coord[1] << " " << coord[2] << endl;
	dx = coord[0] - x0;
	dy = coord[1] - y0;

	cout << "dx: " << dx << ", dy: " << dy << endl;

	/*
	// extract sea bed coordinates
	ZB = new double[nx * ny];
	double pNew[3];
	for (int j = 0; j < dims[1]; j++){
		for (int i = 0; i < dims[0]; i++){
			//cout << i << ", " << j << endl;
			dataset1->GetPoint(i, j, 0, pNew);
			ZB[i * ny + j] = pNew[2];
			cout << pNew[2] << endl;
		}
	}
	*/

	// Calculate beta for all points i LSgrid
	beta = new double[nx * ny * nl];
	double z, welev, seabed, pNew[3];
	for (int k = 0; k < nl; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				//cout << i << ", " << j << endl;
				dataset1->GetPoint(i, j, k, pNew);
				z = pNew[2];
				dataset1->GetPoint(i, j, 0, pNew);
				seabed = pNew[2];
				dataset1->GetPoint(i, j, nl-1, pNew);
				welev = pNew[2];
				beta[k * ny * nx + j * nx + i] = 1. + z2s(z, welev, -seabed);
			}
		}
	}
	

	// Identify field data store in vtk files and see if velocity field can be located. This may either be stored as pointdata or celldata.
	vtkIdType numberOfPointArrays =  dataset1->GetPointData()->GetNumberOfArrays();
	vtkIdType numberOfCellArrays = dataset1->GetCellData()->GetNumberOfArrays();

	if (numberOfPointArrays > 0) {
		// loop over point arrays and identify the index of the velocity array;
		for (vtkIdType i = 0; i < numberOfPointArrays; i++)
		{
			if (strcmp(dataset1->GetPointData()->GetArrayName(i), Uname.c_str()) == 0) {
				Uindex = i;
			}
		}
	}
	if (Uindex == -1) {
		cout << "No point field named " << Uname << " stored in vtk file. Searching for cell array..." << endl;
		for (vtkIdType i = 0; i < numberOfCellArrays; i++)
		{
			if (strcmp(dataset1->GetCellData()->GetArrayName(i), Uname.c_str()) == 0) {
				Uindex = i;
			}
		}
		if (Uindex == -1) {
			cout << "Could not locate specified velocity field " << Uname << " in provided vtk files. Check your waveinput file." << endl;
			exit(-1);
		}
		else {
			cout << "Cell array located with field name " << Uname << ". Cell to point data interpolation switched on." << endl;
			cell2Pointdata = true;
		}
	}

	/*
	// get cell data.
	vtkDataArray* test = dataset1->GetCellData()->GetArray(Uindex);
	vtkIdType numcells = test->GetNumberOfTuples();
	cout <<"numtuples: "<< test->GetNumberOfTuples() << endl;
	
	for (vtkIdType i = 0; i < 10; i++) {
		double* pNe = test->GetTuple3(i);
		cout << pNe[0] << "," << pNe[1] << "," << pNe[2] << "," << endl;
	}
	*/

	if (cell2Pointdata) {
		c2p1 = vtkCellDataToPointData::New();
		// todo: here it is possible to pass only a single celldata field to save time.
		c2p1->PassCellDataOn();
		c2p1->SetInputData(dataset1);
		c2p1->Update();
		U1 = c2p1->GetOutput()->GetPointData()->GetArray(Uindex);
	}
	else {
		U1 = dataset1->GetCellData()->GetArray(Uindex);
	}

	/*
	vtkIdType numpoints = U0->GetNumberOfTuples();
	cout << "numtuples_points: " << U0->GetNumberOfTuples() << endl;
	for (vtkIdType i = 0; i < 10; i++) {
		double* pNe = U0->GetTuple3(i);
		cout << pNe[0] << "," << pNe[1] << "," << pNe[2] << "," << endl;
	}
	*/
	cout << "init complete." << endl;
}



// Function moving data from dataset2 to dataset1, and loads new step into dataset2
void VTKreader::loadNext(string path, const char* fname) {

	vtkDataArray* data;
	// Load file 1
	path.append(fname);

	reader0 = reader1;
	reader0->Update();
	dataset0 = reader0->GetOutput();
	U0 = U1;

	reader1 = vtkXMLStructuredGridReader::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();
	dataset1 = reader1->GetOutput();

	t0 = t1;

	vtkDoubleArray* doubledata = vtkDoubleArray::SafeDownCast(dataset1->GetFieldData()->GetArray("TimeValue"));
	t1 = doubledata->GetValue(0);
	
	cout << "Time interval: " << t0 << " to " << t1 << " sec" << endl;
	dt = t1 - t0;


	if (cell2Pointdata) {
		c2p1 = vtkCellDataToPointData::New();
		c2p1->PassCellDataOn();
		c2p1->SetInputData(dataset1);
		c2p1->Update();
		U1 = c2p1->GetOutput()->GetPointData()->GetArray(Uindex);
	}
	else {
		U1 = dataset1->GetCellData()->GetArray(Uindex);
	}
	loadcount++;
}

template<class T>T* BinarySearch(//<======FIND THE POINTER c | *c <= k < *(c+1)
	T* a, T* b,//<--ARRAY START & END POINTS (ARRAY MUST BE SORTED IN INC. ORDER)
	T k) {//<---------------------------------------------------------------KEY
	if (k < *a)return a - 1;//..........note that a-1 may point to an invalid address
	for (T* c; k < *--b; k > * c ? a = c : b = c + 1)c = a + (b - a) / 2;/*->*/return b;
}

template<class T>T LinInterp(//<===========================LINEAR INTERPOLATOR
	const T* X,//<--------------BRACKETING X VALUES (*X AND X[1] MUST BE VALID)
	const T* Y,//<--------------BRACKETING Y VALUES (*Y AND Y[1] MUST BE VALID)
	T x) {//<--------------VALUE TO INTERPOLATE AT (TYPICALLY, *X <= x <= x[1])
	return*Y + (Y[1] - *Y) * (x - *X) / (X[1] - *X);
}

double VTKreader::stretchInterpLocatorZ(double s, int* iptr, int nxp, int nyp){

	// get z values for given point
	double selevs[nl], laynum[nl];
	for (int k = 0; k < nl; k++) {
		selevs[k] = beta[k* ny * nx + nyp * nx + nxp];
		//cout << "h: " << selevs[k] << endl;
		laynum[k] = double(k);
	}
	// todo: perhaps replace this with a interpolation function from boost at some point
	*iptr = BinarySearch(selevs, selevs + nl, s) - selevs;
	double y = LinInterp(selevs + *iptr, laynum + *iptr, s);

	
	//cout << "s: " << s << ", i: " << *iptr << ", y: " << y << endl;
	return y - double(*iptr);

}


template <typename T>
T clip(const T & n, const T & lower, const T & upper) {
	return std::max(lower, std::min(n, upper));
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double* VTKreader::trilinear_interpolation(double tpt, double xpt, double ypt, double zpt) {

	static double res[5]; // array to store results, where 0 = eta, 1=seabed, 2=u, 3=v, 4=w

	float nxp_temp, nyp_temp;
	double xd = std::modf(clip((xpt - bounds[0]) / dx, 0., nx - 1.), &nxp_temp);
	double yd = std::modf(clip((ypt - bounds[2]) / dy, 0., ny - 1.), &nyp_temp);

	int nxp = int(nxp_temp);
	int nyp = int(nyp_temp);

	double pNew[3];

	// get values of eta
	dataset0->GetPoint(nxp, nyp, nl - 1, pNew);
	double C00 = pNew[2];
	dataset0->GetPoint(nxp, clip(nyp + 1, 0, ny - 1), nl - 1, pNew);
	double C01 = pNew[2];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), nyp, nl - 1, pNew);
	double C10 = pNew[2];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), clip(nyp + 1, 0, ny - 1), nl - 1, pNew);
	double C11 = pNew[2];

	dataset1->GetPoint(nxp, nyp, nl - 1, pNew);
	double D00 = pNew[2];
	dataset1->GetPoint(nxp, clip(nyp + 1, 0, ny - 1), nl - 1, pNew);
	double D01 = pNew[2];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), nyp, nl - 1, pNew);
	double D10 = pNew[2];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), clip(nyp + 1, 0, ny - 1), nl - 1, pNew);
	double D11 = pNew[2];

	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double C0 = C00 * (1. - xd) + C10 * xd;
	double C1 = C01 * (1. - xd) + C11 * xd;
	double D0 = D00 * (1. - xd) + D10 * xd;
	double D1 = D01 * (1. - xd) + D11 * xd;

	double wave_elev0 = C0 * (1. - yd) + C1 * yd;
	double wave_elev1 = D0 * (1. - yd) + D1 * yd;

	// set seabed elevation for point
	res[0] = wave_elev0 * (1. - td) + wave_elev1 * td;

	// get values of seabed
	dataset0->GetPoint(nxp, nyp, 0, pNew);
	C00 = pNew[2];
	dataset0->GetPoint(nxp, clip(nyp + 1, 0, ny - 1), 0, pNew);
	C01 = pNew[2];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), nyp, 0, pNew);
	C10 = pNew[2];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), clip(nyp + 1, 0, ny - 1), 0, pNew);
	C11 = pNew[2];

	dataset1->GetPoint(nxp, nyp, 0, pNew);
	D00 = pNew[2];
	dataset1->GetPoint(nxp, clip(nyp + 1, 0, ny - 1), 0, pNew);
	D01 = pNew[2];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), nyp, 0, pNew);
	D10 = pNew[2];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), clip(nyp + 1, 0, ny - 1), 0, pNew);
	D11 = pNew[2];

	C0 = C00 * (1. - xd) + C10 * xd;
	C1 = C01 * (1. - xd) + C11 * xd;
	D0 = D00 * (1. - xd) + D10 * xd;
	D1 = D01 * (1. - xd) + D11 * xd;

	double zb0 = C0 * (1. - yd) + C1 * yd;
	double zb1 = D0 * (1. - yd) + D1 * yd;

	// set seabed elevation for point
	res[1] =  zb0 * (1. - td) + zb1 * td;


	//std::cout << "dy: " << dy << ", nyp: " << nyp << ", ypt: " << ypt << std::endl;

	double spt0 = std::max(z2s(std::min(zpt, wave_elev0), wave_elev0, -zb0), -1.);
	double spt1 = std::max(z2s(std::min(zpt, wave_elev1), wave_elev1, -zb0), -1.);
	

	int nsp0, nsp1;

	
	double sd0 = stretchInterpLocatorZ(spt0 + 1., &nsp0, nxp, nyp);
	double sd1 = stretchInterpLocatorZ(spt1 + 1., &nsp1, nxp, nyp);

	int temp = clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1);
	
	//cout << "nxp: " << nxp << ", nyp: " << nyp << "nsp0: " << nsp0 << "sum: " << temp << endl;

	//exit(0);
	// Trilinear interpolation.
	double* C000 = U0->GetTuple3(nxp * ny * nl + nyp * nl + nsp0);

	double* C001 = U0->GetTuple3(nxp * ny * nl + nyp * nl + clip(nsp0 + 1, 0, nl - 1));

	double* C010 = U0->GetTuple3(nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp0);

	double* C011 = U0->GetTuple3(nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1));

	double* C100 = U0->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nsp0);

	double* C101 = U0->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clip(nsp0 + 1, 0, nl - 1));

	double* C110 = U0->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp0);

	double* C111 = U0->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1));
	

	double* D000 = U1->GetTuple3(nxp * ny * nl + nyp * nl + nsp1);

	double* D001 = U1->GetTuple3(nxp * ny * nl + nyp * nl + clip(nsp1 + 1, 0, nl - 1));

	double* D010 = U1->GetTuple3(nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp1);

	double* D011 = U1->GetTuple3(nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp1 + 1, 0, nl - 1));

	double* D100 = U1->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nsp1);

	double* D101 = U1->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clip(nsp1 + 1, 0, nl - 1));

	double* D110 = U1->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp1);

	double* D111 = U1->GetTuple3(clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp1 + 1, 0, nl - 1));

	//double sd0 = nsp0 - floor(nsp0);
	//double sd1 = nsp1 - floor(nsp1);

	double CC00[3], CC01[3], CC10[3], CC11[3];
	CC00[0] = C000[0] * (1. - xd) + C100[0] * xd;
	CC00[1] = C000[1] * (1. - xd) + C100[1] * xd;
	CC00[2] = C000[2] * (1. - xd) + C100[2] * xd;

	CC01[0] = C001[0] * (1. - xd) + C101[0] * xd;
	CC01[1] = C001[1] * (1. - xd) + C101[1] * xd;
	CC01[2] = C001[2] * (1. - xd) + C101[2] * xd;

	CC10[0] = C010[0] * (1. - xd) + C110[0] * xd;
	CC10[1] = C010[1] * (1. - xd) + C110[1] * xd;
	CC10[2] = C010[2] * (1. - xd) + C110[2] * xd;

	CC11[0] = C011[0] * (1. - xd) + C111[0] * xd;
	CC11[1] = C011[1] * (1. - xd) + C111[1] * xd;
	CC11[2] = C011[2] * (1. - xd) + C111[2] * xd;

	//std::cout << int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0)) << ", " << ceil(nxp) << ", " << ceil(nyp) << ", " << ceil(nsp0) << std::endl;

	//std::cout << C010 << ", " << C110 << ", " << C011 << ", " << C111 << std::endl;

	double DD00[3], DD01[3], DD10[3], DD11[3];
	DD00[0] = D000[0] * (1. - xd) + D100[0] * xd;
	DD00[1] = D000[1] * (1. - xd) + D100[1] * xd;
	DD00[2] = D000[2] * (1. - xd) + D100[2] * xd;
	DD01[0] = D001[0] * (1. - xd) + D101[0] * xd;
	DD01[1] = D001[1] * (1. - xd) + D101[1] * xd;
	DD01[2] = D001[2] * (1. - xd) + D101[2] * xd;
	DD10[0] = D010[0] * (1. - xd) + D110[0] * xd;
	DD10[1] = D010[1] * (1. - xd) + D110[1] * xd;
	DD10[2] = D010[2] * (1. - xd) + D110[2] * xd;
	DD11[0] = D011[0] * (1. - xd) + D111[0] * xd;
	DD11[1] = D011[1] * (1. - xd) + D111[1] * xd;
	DD11[2] = D011[2] * (1. - xd) + D111[2] * xd;

	double CC0[3], CC1[3], DD0[3], DD1[3];
	CC0[0] = CC00[0] * (1. - yd) + CC10[0] * yd;
	CC0[1] = CC00[1] * (1. - yd) + CC10[1] * yd;
	CC0[2] = CC00[2] * (1. - yd) + CC10[2] * yd;
	CC1[0] = CC01[0] * (1. - yd) + CC11[0] * yd;
	CC1[1] = CC01[1] * (1. - yd) + CC11[1] * yd;
	CC1[2] = CC01[2] * (1. - yd) + CC11[2] * yd;
	DD0[0] = DD00[0] * (1. - yd) + DD10[0] * yd;
	DD0[1] = DD00[1] * (1. - yd) + DD10[1] * yd;
	DD0[2] = DD00[2] * (1. - yd) + DD10[2] * yd;
	DD1[0] = DD01[0] * (1. - yd) + DD11[0] * yd;
	DD1[1] = DD01[1] * (1. - yd) + DD11[1] * yd;
	DD1[2] = DD01[2] * (1. - yd) + DD11[2] * yd;

	
	res[2] = (CC0[0] * (1. - sd0) + CC1[0] * sd0) * (1. - td) + (DD0[0] * (1. - sd1) + DD1[0] * sd1) * td; // u
	res[3] = (CC0[1] * (1. - sd0) + CC1[1] * sd0) * (1. - td) + (DD0[1] * (1. - sd1) + DD1[1] * sd1) * td; // v
	res[4] = (CC0[2] * (1. - sd0) + CC1[2] * sd0) * (1. - td) + (DD0[2] * (1. - sd1) + DD1[2] * sd1) * td; // w

	return res;
}

bool VTKreader::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if (tpt > t1) {
		std::cout << "t0: " << t0 << ", t1: " << (t1) << ", tpt: " << tpt << std::endl;
		return false;

	}
	return true;
}

//----------------------------------------------------------------------------------------------------------------------------------------


// VTK output


//----------------------------------------------------------------------------------------------------------------------------------------


// when called, writes stored kinematics to file
void VTKreader::write_vtk(bool endtime) {
	
	char buffer[256];
	if (endtime) {
		sprintf(buffer, "%05d", loadcount);
	}
	else {
		sprintf(buffer, "%05d", loadcount - 1);
	}


	std::string vtk_directory_path = "./";
	
	if (dirExists(vtk_directory_path.c_str()) == 0) {
		std::cout << "WARNING: Specified directory for storage of VTK files does not exist. Directory will be created at the following path:  " << vtk_directory_path << std::endl;
		createDirectory(vtk_directory_path);
	}

	std::string str(buffer);
	std::string fpath = (vtk_directory_path + vtk_prefix + buffer + ".vtu");
	std::cout << fpath << std::endl;
	FILE* fp = fopen(fpath.c_str(), "w");
	if (endtime)
		export_vtu(fp, true);
	else
		export_vtu(fp, false);
	fclose(fp);

	std::cout << "wrote kinematics to: " << fpath << std::endl;
}

/* exports read vts data to vtu (For QA purposes) */
void VTKreader::export_vtu(FILE* fp, bool last)
{
	// write header
	fputs("<?xml version=\"1.0\"?>\n"
		"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
	fputs("\t <UnstructuredGrid>\n", fp);
	fprintf(fp, "\t\t <FieldData> \n");
	if (last) {
		fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", vtk_timelabel.c_str(), t0 + dt, t0 + dt);
		fprintf(fp, "\t\t\t %.3f \n", t0 + dt);
	}
	else {
		fprintf(fp, "\t\t\t <DataArray type = \"Float64\" Name = \"%s\" NumberOfTuples = \"1\" format = \"ascii\" RangeMin = \"%.3f\" RangeMax = \"%.3f\"> \n", vtk_timelabel.c_str(), t0, t0);
		fprintf(fp, "\t\t\t %.3f \n", t0);
	}
	fprintf(fp, "\t\t\t </DataArray > \n");
	fprintf(fp, "\t\t </FieldData> \n");

	fprintf(fp, "\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nx * ny * nl, std::max((nx - 1), 1) * std::max((ny - 1), 1) * std::max((nl - 1), 1));

	// Loop over velocity data and store kinematics in cell vector stucture
	fputs("\t\t\t <PointData Scalars=\"scalars\">\n", fp);

	fprintf(fp, "\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">\n");
	if (last) {
		for (int m = 0; m < nl; m++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {	
					double* C1 = U1->GetTuple3(m * ny * nx + j * nx + i);
					fprintf(fp, "%g %g %g\n", C1[0], C1[1], C1[2]);
				}
			}
		}
	}
	else {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int m = 0; m < nl; m++) {
					double* C0 = U0->GetTuple3(m * ny * nx + j * nx + i);
					fprintf(fp, "%g %g %g\n", C0[0], C0[1], C0[2]);
				}
			}
		}
	}
	fputs("\t\t\t\t </DataArray>\n", fp);

	fputs("\t\t\t </PointData>\n", fp);

	fputs("\t\t\t <Points>\n", fp);
	fputs("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
	double pNew[3];
	if (last) {
		for (int i = 0; i < nx; i++) {
			double xpt = bounds[0] + dx * i;
			for (int j = 0; j < ny; j++) {
				double ypt = bounds[2] + dy * j;
				// get values of eta
				dataset1->GetPoint(i, j, nl - 1, pNew);
				double eta1_temp = pNew[2];
				dataset1->GetPoint(i, j, 0, pNew);
				double seabed = pNew[2];
				for (int m = 0; m < nl; m++) {
					double zpt1 = s2z(beta[i * ny * nl + j * nl + m] - 1., eta1_temp, -seabed);
					fprintf(fp, "%12.4f %12.4f %12.4f\n", xpt, ypt, zpt1);
				}
			}
		}
	}
	else {
		for (int i = 0; i < nx; i++) {
			double xpt = bounds[0] + dx * i;
			for (int j = 0; j < ny; j++) {
				double ypt = bounds[2] + dy * j;
				// get values of eta
				dataset0->GetPoint(i, j, nl - 1, pNew);
				double eta0_temp = pNew[2];
				//ypt = pNew[1];
				//xpt = pNew[0];
				dataset0->GetPoint(i, j, 0, pNew);
				double seabed = pNew[2];
				for (int m = 0; m < nl; m++) {
					double zpt0 = s2z(beta[m * ny * nx + j * nx + i] - 1., eta0_temp, -seabed);
					fprintf(fp, "%12.4f %12.4f %12.4f\n", xpt, ypt, zpt0);
				}
			}
		}
	}
	fputs("\t\t\t\t </DataArray>\n", fp);
	fputs("\t\t\t </Points>\n", fp);

	fputs("\t\t\t <Cells>\n", fp);
	fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);

	if (nx > 1 && ny > 1 && nl > 1) {
		for (int i = 0; i < (nx - 1); i++) {
			for (int j = 0; j < (ny - 1); j++) {
				for (int m = 0; m < (nl - 1); m++) {
					int ape1 = nl * ny * i + nl * j + m;
					int ape2 = nl * ny * (i + 1) + nl * j + m;
					int ape3 = nl * ny * (i + 1) + nl * (j + 1) + m;
					int ape4 = nl * ny * i + nl * (j + 1) + m;
					int ape5 = nl * ny * i + nl * j + (m + 1);
					int ape6 = nl * ny * (i + 1) + nl * j + (m + 1);
					int ape7 = nl * ny * (i + 1) + nl * (j + 1) + (m + 1);
					int ape8 = nl * ny * i + nl * (j + 1) + (m + 1);
					fprintf(fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);
				}
			}
		}


		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int i = 1; i < ((nx - 1) * (ny - 1) * (nl - 1) + 1); i++) {
			fprintf(fp, "%d \n", i * 8);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int i = 1; i < ((nx - 1) * (ny - 1) * (nl - 1) + 1); i++) {
			fputs("12 \n", fp);
		}
	}
	// only single dimension i y direction.
	else if (nx > 1 && ny == 1 && nl > 1) {
		for (int i = 0; i < (nx - 1); i++) {
			for (int m = 0; m < (nl - 1); m++) {
				int ape1 = nl * i + m;
				int ape2 = nl * (i + 1) + m;
				int ape3 = nl * (i + 1) + (m + 1);
				int ape4 = nl * i + (m + 1);
				fprintf(fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4);
			}

		}

		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int i = 1; i < ((nx - 1) * (nl - 1) + 1); i++) {
			fprintf(fp, "%d \n", i * 4);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int i = 1; i < ((nx - 1) * (nl - 1) + 1); i++) {
			fputs("9 \n", fp);
		}

	}
	// only single dimension i x direction.
	else if (nx == 1 && ny > 1 && nl > 1) {
		for (int j = 0; j < (ny - 1); j++) {
			for (int m = 0; m < (nl - 1); m++) {
				int ape1 = nl * j + m;
				int ape2 = nl * (j + 1) + m;
				int ape3 = nl * (j + 1) + (m + 1);
				int ape4 = nl * j + (m + 1);
				fprintf(fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4);
			}

		}

		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int j = 1; j < ((ny - 1) * (nl - 1) + 1); j++) {
			fprintf(fp, "%d \n", j * 4);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int j = 1; j < ((ny - 1) * (nl - 1) + 1); j++) {
			fputs("9 \n", fp);
		}
	}

	// only single dimension i z direction (lagrangian).
	else if (nx > 1 && ny > 1 && nl == 1) {
		for (int i = 0; i < (nx - 1); i++) {
			for (int j = 0; j < (ny - 1); j++) {
				int ape1 = ny * i + j;
				int ape2 = ny * (i + 1) + j;
				int ape3 = ny * (i + 1) + (j + 1);
				int ape4 = ny * i + (j + 1);
				fprintf(fp, "%u %u %u %u\n", ape1, ape2, ape3, ape4);
			}

		}

		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

		for (int i = 1; i < ((nx - 1) * (ny - 1) + 1); i++) {
			fprintf(fp, "%d \n", i * 4);
		}
		fputs("\t\t\t\t </DataArray>\n", fp);
		fputs("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
		for (int i = 1; i < ((nx - 1) * (ny - 1) + 1); i++) {
			fputs("9 \n", fp);
		}

	}

	fputs("\t\t\t\t </DataArray>\n", fp);
	fputs("\t\t\t </Cells>\n", fp);
	fputs("\t\t </Piece>\n", fp);
	fputs("\t </UnstructuredGrid>\n", fp);
	fputs("</VTKFile>\n", fp);
	fflush(fp);
}

