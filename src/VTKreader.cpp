#include "VTKreader.h"
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
		trilinear_interpolation(0.1, 130., 0.1, 0.2);
		//loadNext(vtkfilepath, filevec->at(1).c_str());
		filecount = 1;
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

	reader0 = vtkXMLStructuredGridReader::New();
	reader0->SetFileName(path.c_str());
	reader0->Update();
	dataset0 = reader0->GetOutput();

	reader1 = vtkXMLStructuredGridReader::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();
	dataset1 = reader1->GetOutput();
	

	// Get time
	vtkDataArray* data;
	data = dataset0->GetFieldData()->GetArray("TimeValue");
	doubledata = vtkDoubleArray::SafeDownCast(data);
	t0 = doubledata->GetValue(0);
	data = dataset1->GetFieldData()->GetArray("TimeValue");
	doubledata = vtkDoubleArray::SafeDownCast(data);
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
	dataset1->GetPoint(1, 1, 0, coord);
	dx = coord[0] - x0;
	dy = coord[1] - y0;


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

	// Calculate beta for all points i LSgrid
	beta = new double[nx * ny * nl];
	double z, welev, seabed;
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

				beta[i * ny * nl + j * nl + k] = 1. + z2s(z, welev, -seabed);
			}
		}
	}
	/*
	// extract surface coordinates
	ETA1 = new double[nx * ny];
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			dataset1->GetPoint(i, j, nl-1, pNew);
			ETA1[i * ny + j] = pNew[2];
		}
	}
	*/

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


	// get cell data.
	vtkDataArray* test = dataset1->GetCellData()->GetArray(Uindex);
	vtkIdType numcells = test->GetNumberOfTuples();
	cout <<"numtuples: "<< test->GetNumberOfTuples() << endl;
	
	for (vtkIdType i = 0; i < 10; i++) {
		double* pNe = test->GetTuple3(i);
		cout << pNe[0] << "," << pNe[1] << "," << pNe[2] << "," << endl;
	}

	// cell 2 pointdata
	vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();
	// todo: here it is possible to pass only a single celldata field to save time.
	c2p->PassCellDataOn();
	c2p->SetInputData(dataset1);
	c2p->Update();
	vtkDataArray* test2 = c2p->GetOutput()->GetPointData()->GetArray(Uindex);

	vtkIdType numpoints = test2->GetNumberOfTuples();
	cout << "numtuples_points: " << test2->GetNumberOfTuples() << endl;
	for (vtkIdType i = 0; i < 10; i++) {
		double* pNe = test2->GetTuple3(i);
		cout << pNe[0] << "," << pNe[1] << "," << pNe[2] << "," << endl;
	}
	
	
	/* steps to perform:
	 1. check if vector field velocity exists as point data if not
	 2. check if vector field velocity exists as cell data, is so, convert to pointdata switch on.
	 3. extract point data
	 3b. extract ncellx, y, z
	 3c. extract bottom coordinates
	 3d. extract surface coordinates
	 3e.perform inverse sigma transform and calculate stretch factor for each point
	 4. download velocity data

	*/

}

/*
double VTKreader::locatorZ() {

	/*
	1. find locator horizontally (x and y)
	2. find
	

}
*/

template <typename T>
T clip(const T & n, const T & lower, const T & upper) {
	return std::max(lower, std::min(n, upper));
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double* VTKreader::trilinear_interpolation(double tpt, double xpt, double ypt, double zpt) {

	static double res[4]; // array to store results, where 0 = eta, 1=u, 2=v, 3=w

	float nxp_temp, nyp_temp;
	double xd = std::modf(clip((xpt - bounds[0]) / dx, 0., nx - 1.), &nxp_temp);
	double yd = std::modf(clip((ypt - bounds[2]) / dy, 0., ny - 1.), &nyp_temp);

	int nxp = int(nxp_temp);
	int nyp = int(nyp_temp);

	double pNew[3];

	dataset0->GetPoint(nxp, nyp, nl - 1, pNew);
	double C00 = pNew[2];
	dataset0->GetPoint(nxp, clip(nyp + 1, 0, ny - 1), nl - 1, pNew);
	double C01 = pNew[2];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), nyp, nl - 1, pNew);
	double C10 = pNew[2];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), clip(nyp + 1, 0, ny - 1), nl - 1, pNew);
	double C11 = pNew[2];

	//cout << "nxp: " << nxp << ", nyp:" << nyp << endl;
	//cout << " C00: " << C00 << " C01: " << C01 << " C10: " << C10 << " C11: " << C00 << endl;

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

	// set wave elevation for point
	res[0] =  (wave_elev0 * (1. - td) + wave_elev1 * td;

	//cout << "wave elev0: " << wave_elev0 << ", wave elev1: " << wave_elev1 << endl;
	
	
	
	


	//std::cout << "dy: " << dy << ", nyp: " << nyp << ", ypt: " << ypt << std::endl;

	double spt0 = std::max(z2s(std::min(zpt - swl, wave_elev0 - swl), wave_elev0 - swl, water_depth), -1.);
	double nsp0a = locator
	double nsp0a = (spt0 + 1.) / ds;
	double spt1 = std::max(z2s(std::min(zpt - swl, wave_elev1 - swl), wave_elev1 - swl, water_depth), -1.);
	double nsp1a = (spt1 + 1.) / ds;

	float nsp0_temp, nsp1_temp;
	double sd0 = std::modf((spt0 + 1.) / ds, &nsp0_temp);
	double sd1 = std::modf((spt1 + 1.) / ds, &nsp1_temp);
	int nsp0 = int(nsp0_temp);
	int nsp1 = int(nsp1_temp);

	//std::cout << floor(nsp0a) << " " << nsp0 << "spt0: " << spt0 << ", wavelev0: " << wave_elev0 << std::endl;

	//std::cout << "nsp0:" << nsp0 << ", nsp1: " << nsp1 << std::endl;

	//exit(0);
	// Trilinear interpolation.
	double C000 = VAR0[nxp * ny * nl + nyp * nl + nsp0];
	double C001 = VAR0[nxp * ny * nl + nyp * nl + clip(nsp0 + 1, 0, nl - 1)];
	double C010 = VAR0[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp0];
	double C011 = VAR0[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1)];
	double C100 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nsp0];
	double C101 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clip(nsp0 + 1, 0, nl - 1)];
	double C110 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp0];//
	double C111 = VAR0[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1)];//
	double D000 = VAR1[nxp * ny * nl + nyp * nl + nsp1];
	double D001 = VAR1[nxp * ny * nl + nyp * nl + clip(nsp1 + 1, 0, nl - 1)];
	double D010 = VAR1[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp1];
	double D011 = VAR1[nxp * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp1 + 1, 0, nl - 1)];
	double D100 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + nsp1];
	double D101 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + nyp * nl + clip(nsp1 + 1, 0, nl - 1)];
	double D110 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + nsp1];
	double D111 = VAR1[clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp1 + 1, 0, nl - 1)];

	//double sd0 = nsp0 - floor(nsp0);
	//double sd1 = nsp1 - floor(nsp1);

	C00 = C000 * (1. - xd) + C100 * xd;
	C01 = C001 * (1. - xd) + C101 * xd;
	C10 = C010 * (1. - xd) + C110 * xd;
	C11 = C011 * (1. - xd) + C111 * xd;
	//std::cout << int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0)) << ", " << ceil(nxp) << ", " << ceil(nyp) << ", " << ceil(nsp0) << std::endl;

	//std::cout << C010 << ", " << C110 << ", " << C011 << ", " << C111 << std::endl;

	D00 = D000 * (1. - xd) + D100 * xd;
	D01 = D001 * (1. - xd) + D101 * xd;
	D10 = D010 * (1. - xd) + D110 * xd;
	D11 = D011 * (1. - xd) + D111 * xd;

	C0 = C00 * (1. - yd) + C10 * yd;
	C1 = C01 * (1. - yd) + C11 * yd;
	D0 = D00 * (1. - yd) + D10 * yd;
	D1 = D01 * (1. - yd) + D11 * yd;

	//return (C0 * (1. - sd0) + C1 * sd0) * (1. - td) + (D0 * (1. - sd1) + D1 * sd1) * td;

	

	return res;
}


/*
// Function moving data from dataset2 to dataset1, and loads new step into dataset2
void VTKreader::loadNext(string path, const char* fname) {

	vtkDataArray* data;
	// Load file 1
	path.append(fname);
	//cout << path << endl;

	reader0 = reader1;
	reader0->Update();

	dataset0 = reader0->GetOutput(); // move dataset1 to dataset0	
	
					 //dataset0 = dataset1;
	//reader1 = vtkXMLStructuredGridReader::New();
	reader1 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();

	dataset1 = reader1->GetOutput(); // load new dataset into dataset1
									 //cout << "loaded successfully!\n";

									 // load time variable from the datasets
									 // Extract Time from dataset of file 1
	data = dataset0->GetFieldData()->GetArray("TIME");
	doubledata = vtkDoubleArray::SafeDownCast(data);
	time0 = doubledata->GetValue(0);
	data = dataset1->GetFieldData()->GetArray("TIME");
	doubledata = vtkDoubleArray::SafeDownCast(data);
	time1 = doubledata->GetValue(0);
	cout << "Time interval: " << time0 << " to " << time1 << " sec" << endl;
	//data->Delete();
	//cstring vfraqstr(chF.c_str());
	//cstring velostr(chV.c_str());
	// Check that the later timestep file contains the field variables needed


	vtkIdType numberOfPointArrays = dataset1->GetPointData()->GetNumberOfArrays();
	vfraq_field_located = 0;
	velo_field_located = 0;
	for (vtkIdType i = 0; i < numberOfPointArrays; i++)
	{

		if (strcmp(dataset1->GetPointData()->GetArrayName(i), Uname.c_str()) == 0) {
			velo_field_located = 1;
		}
	}
}
*/
