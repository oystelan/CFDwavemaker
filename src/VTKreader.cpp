#include "VTKreader.h"
#include "Utils.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

// ---------------------------------------------------------------------------
// Piecewise-parabolic (PPR) reconstruction of a cell-centred vertical profile.
// Ported from the remap_central kernel of ~/progs/basilisk/src/layered/remapc.h:
// a parabola is fitted per layer honouring the layer integral (cell average) and
// the two interface values, with monotonic limiting -> conservative and monotone.
//
//   zi[0..ncell]     : interface z positions (increasing), zi[k] lower face of cell k
//   zc[0..ncell-1]   : cell-centre z positions
//   val[0..ncell-1]  : cell-averaged values (velocity component)
//   h[0..ncell-1]    : layer thicknesses (= zi[k+1]-zi[k])
//   zq               : query height,   method: 0 = PPR (default), 1 = linear between centres
// ---------------------------------------------------------------------------
static double ppr_reconstruct(const double* zi, const double* zc, const double* val,
                              const double* h, int ncell, double zq, int method)
{
	if (ncell <= 0) return 0.0;
	if (ncell == 1) return val[0];

	// clamp query to the column
	if (zq < zi[0]) zq = zi[0];
	if (zq > zi[ncell]) zq = zi[ncell];

	if (method == 1) {
		// linear between cell centres (fallback)
		if (zq <= zc[0]) return val[0];
		if (zq >= zc[ncell - 1]) return val[ncell - 1];
		int kk = 0;
		while (kk < ncell - 2 && zq > zc[kk + 1]) kk++;
		double wgt = (zq - zc[kk]) / (zc[kk + 1] - zc[kk]);
		return val[kk] * (1.0 - wgt) + val[kk + 1] * wgt;
	}

	// locate layer k containing zq
	int k = 0;
	while (k < ncell - 1 && zq > zi[k + 1]) k++;

	// interface values at zi[k] (sl) and zi[k+1] (sr), h-weighted linear between centres
	auto iface = [&](int m) -> double {
		if (m <= 0)        return 1.5 * val[0] - 0.5 * val[1];
		if (m >= ncell)    return 1.5 * val[ncell - 1] - 0.5 * val[ncell - 2];
		double t = (zi[m] - zc[m - 1]) / (zc[m] - zc[m - 1]);
		return val[m - 1] + t * (val[m] - val[m - 1]);
	};
	double sl = iface(k), sr = iface(k + 1);
	double fk = val[k];

	// monotonic limiting of sl, sr (remap_central), only where neighbours exist
	if (k > 0 && k < ncell - 1) {
		if ((val[k + 1] - fk) * (fk - val[k - 1]) < 0.0) {
			sl = fk; sr = fk;
		} else {
			double dx = h[k];
			double sigma_left  = 2.0 * (fk - val[k - 1]) / dx;
			double sigma_center= 2.0 * (val[k + 1] - val[k - 1]) / (zi[k + 2] - zi[k - 1] + dx);
			double sigma_right = 2.0 * (val[k + 1] - fk) / dx;
			double sabs = std::min(std::fabs(sigma_center), std::min(std::fabs(sigma_left), std::fabs(sigma_right)));
			double sigma = (sigma_left * sigma_right <= 0.0 || sigma_center * sigma_left <= 0.0) ? 0.0
			               : (sigma_left < 0.0 ? -sabs : sabs);
			if ((sl - val[k - 1]) * (fk - sl) < 0.0) sl = fk - 0.5 * dx * sigma;
			if ((sr - fk) * (val[k + 1] - sr) < 0.0) sr = fk + 0.5 * dx * sigma;
		}
	}

	// parabola honouring cell average fk:  P(z) = C0 + C1 zeta + C2 zeta^2,  zeta in [0,1]
	double C0 = sl, C1 = 6.0 * fk - 2.0 * sr - 4.0 * sl, C2 = 3.0 * (sr + sl - 2.0 * fk);
	// clip extremum to avoid over/undershoot
	if (C1 * C2 < 0.0 && C1 / C2 > -2.0) {
		if (C1 / C2 > -1.0) sr = 3.0 * fk - 2.0 * sl;
		else                sl = 3.0 * fk - 2.0 * sr;
		C0 = sl; C1 = 6.0 * fk - 2.0 * sr - 4.0 * sl; C2 = 3.0 * (sr + sl - 2.0 * fk);
	}
	double zeta = (zq - zi[k]) / h[k];
	return C0 + C1 * zeta + C2 * zeta * zeta;
}

// Transforms normal z axis to streched sigma coordinates 
// s defined between -1 (seabed) and 0 (free surface)
double VTKreader::z2s(double z, double wave_elev, double depth) {
	
	//cout << "z2s: " << z << " " << wave_elev << " " << depth << endl;
	//return (z - wave_elev) / max((depth + wave_elev),0.00000001);
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
void VTKreader::init(double tpt) {

	filevec = listdir(vtkfilepath.c_str(), vtk_prefix.c_str());
	ostream_iterator<string> out(cout, " ");
	cout << "The following list of files identified in specified folder:" << endl;
	copy(filevec->begin(), filevec->end(), out);
	cout << endl;

	//cout << filevec->size();
	if (filevec->size() >= 2) {
		// get min and max time value store in fileseries

		tmin = getTimeFromVTKFile(vtkfilepath, filevec->at(0).c_str());
		tmax = getTimeFromVTKFile(vtkfilepath, filevec->at(filevec->size()-1).c_str());
		
		cout << "Time range of vts files: " << tmin << " to " << tmax << " seconds." << endl;

		// calculate which file to start from
		double dt_approx = (tmax - tmin) / (filevec->size() - 1);


		loadcount = max(int(floor((tpt + t_start - tmin) / dt_approx))-1, 0);

		cout << "Approximated file startfile no. " << loadcount << endl;

		// load vtufiles
		loadInit(vtkfilepath, filevec->at(loadcount).c_str());
		loadNext(vtkfilepath, filevec->at(loadcount+1).c_str());
		//write_vtk(false);
		cout << "VTU Files found and loaded..." << endl;
		//cout << "u: " << u(0.1, 70., 0.0, 0.0) << endl;
	}
	else {
		cout << "Not enough (or no) VTU files in folder. A minimum of two files are required.\n\
			Perhaps the wrong path is specified? just trying to be helpful..." << endl;
	}

}


void VTKreader::update(double tpt) {

	if (loadcount < (filevec->size() + 1)) {
		int i = 1;
		while (loadcount < (filevec->size() - 1)) {
			loadNext(vtkfilepath, filevec->at(loadcount + 1).c_str());
			if ((tpt + t_start) >= t0 && (tpt + t_start) <= t1) {
				break;
			}
		}
	}
	else {
		cerr << "End of VTK file list reached. (i.e. no more data to load)." << endl;
		exit(-2);
	}
}


// Full kinematics at (t,x,y,z) with a thread-local one-entry cache. Because the
// solver typically asks for VeloX, VeloY, VeloZ and nhPres at the same point in
// succession, caching the last evaluation avoids recomputing the (expensive)
// interpolation once per component.
void VTKreader::get_kinematics(double tpt, double xpt, double ypt, double zpt, double* out) {
	static thread_local double ct = -1e300, cx = 0., cy = 0., cz = 0.;
	static thread_local double cres[6] = { 0.,0.,0.,0.,0.,0. };
	static thread_local bool cvalid = false;

	if (cvalid && tpt == ct && xpt == cx && ypt == cy && zpt == cz) {
		for (int i = 0; i < 6; i++) out[i] = cres[i];
		return;
	}
	if (has_h) {
		multilayer_interpolation(cres, tpt + t_start, xpt, ypt, zpt);
	}
	else if (input2d) {
		bilinear_interpolation_xy(cres, tpt + t_start, xpt, zpt);
		cres[5] = 0.0;
	}
	else {
		trilinear_interpolation(cres, tpt + t_start, xpt, ypt, zpt);
		cres[5] = 0.0;
	}
	ct = tpt; cx = xpt; cy = ypt; cz = zpt; cvalid = true;
	for (int i = 0; i < 6; i++) out[i] = cres[i];
}

double VTKreader::u(double tpt, double xpt, double ypt, double zpt) {
	double r[6]; get_kinematics(tpt, xpt, ypt, zpt, r); return r[2];
}

double VTKreader::v(double tpt, double xpt, double ypt, double zpt) {
	double r[6]; get_kinematics(tpt, xpt, ypt, zpt, r); return r[3];
}

double VTKreader::w(double tpt, double xpt, double ypt, double zpt) {
	double r[6]; get_kinematics(tpt, xpt, ypt, zpt, r); return r[4];
}

double VTKreader::eta(double tpt, double xpt, double ypt) {
		double res[6];
		double* temp;
		if (has_h) { double r[6]; get_kinematics(tpt, xpt, ypt, 0., r); return r[0]; }
        if (input2d){
			 temp = linear_interpolation(res, (tpt + t_start), xpt);
        }
        else{
			 temp = bilinear_interpolation(res, (tpt + t_start), xpt, ypt);
        }
        return temp[0];
}

double VTKreader::seabed(double xpt, double ypt) {
		double res[2];
		double* temp;
        if (input2d){
			 temp = linear_interpolation(res, 0, xpt);
        }
        else{
			 temp = bilinear_interpolation(res, 0, xpt, ypt);
        }
	return temp[1];
}


double VTKreader::getTimeFromVTKFile(string path, const char* fname) {
	path.append(fname);
	vtkXMLStructuredGridReader* tempreader = vtkXMLStructuredGridReader::New();
	tempreader->SetFileName(path.c_str());
	tempreader->Update();
	vtkStructuredGrid* dataset = tempreader->GetOutput();

	// Get time
	vtkDoubleArray* doubledata = vtkDoubleArray::SafeDownCast(dataset->GetFieldData()->GetArray("TimeValue"));
	
	return doubledata->GetValue(0);
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
	
	
	//get extent of grid (ncells in each direction
	int* dims = dataset1->GetDimensions();
	cout << "Dimensions of grid: " << dims[0] << ", " << dims[1] << ", " << dims[2] << ", " << endl;

	// check if data is 2D or 3D
	if (dims[2] == 1) {
		cout << "Two dimensional data read, given in x-y plane" << endl;
		input2d = true;
		zmin = bounds[2];
	}
	else {
		cout << "This is a 3D case." << endl;
		zmin = bounds[4];
	}


	if (input2d) {
		nx = dims[0];
		ny = 1;
		nl = dims[1];
	}
	else {
		nx = dims[0];
		ny = dims[1];
		nl = dims[2];
	}

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

	// Compute height function betah
	if (input2d) {
		// Calculate betah for all points i LSgrid
		betah = new double[nx * nl];
		double z, welev, seabed, pNew[3];
		for (int k = 0; k < nl; k++) {
			for (int i = 0; i < nx; i++) {
				//cout << i << ", " << j << endl;
				dataset1->GetPoint(i, k, 0, pNew);
				z = pNew[1];
				dataset1->GetPoint(i, 0, 0, pNew);
				seabed = pNew[1];
				dataset1->GetPoint(i, nl - 1, 0, pNew);
				welev = pNew[1];
				betah[k * nx  + i] = welev == seabed ? 0. : 1 + z2s(z, welev, -seabed);
				//cout << "seabed:" << seabed << " welev: " << welev << " betah:" << betah[k * nx  + i]<< endl;
			}
		}
	}
	else {
		//cout << "nx, ny, nl: " << nx << ", " << ny << ", " << nl << endl;
		// Calculate betah for all points i LSgrid
		betah = new double[nx * ny * nl];
		double z, welev, seabed, pNew[3];
		for (int k = 0; k < nl; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					//cout << i << ", " << j << ", " << k << endl;
					dataset1->GetPoint(i, j, k, pNew);
					z = pNew[2];
					dataset1->GetPoint(i, j, 0, pNew);
					seabed = pNew[2];
					dataset1->GetPoint(i, j, nl - 1, pNew);
					welev = pNew[2];
					
					betah[k * ny * nx + j * nx + i] = welev == seabed ? 0. : 1. + z2s(z, welev, -seabed);
					//cout << "seabed:" << seabed << " welev: " << welev << " z: " << z << " betah:" << betah[k * ny * nx + j * nx + i] << endl;
				}
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

	// Detect optional multilayer cell fields (h, phi). If h is present the reader
	// switches to the h-based reconstruction and uses raw cell-centred velocity
	// (no vtkCellDataToPointData averaging).
	locate_multilayer_fields();
	if (has_h) {
		cout << "Multilayer cell field 'h' located: eta and velocity reconstructed from layer thicknesses." << endl;
		if (has_phi) cout << "Non-hydrostatic pressure field 'phi' located (defined on lower layer interfaces)." << endl;
		cout << "Vertical velocity reconstruction: " << (vertical_interp == 0 ? "PPR (parabolic)" : "linear") << endl;
	}

	if (has_h) {
		// raw cell velocity already bound in locate_multilayer_fields (Ucell1); U1 unused
	}
	else if (cell2Pointdata) {
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
	// carry the multilayer cell arrays over to the "previous" timestep slot
	Ucell0 = Ucell1;
	Hcell0 = Hcell1;
	Phicell0 = Phicell1;

	reader1 = vtkXMLStructuredGridReader::New();
	reader1->SetFileName(path.c_str());
	reader1->Update();
	dataset1 = reader1->GetOutput();

	t0 = t1;

	vtkDoubleArray* doubledata = vtkDoubleArray::SafeDownCast(dataset1->GetFieldData()->GetArray("TimeValue"));
	t1 = doubledata->GetValue(0);

	cout << "Loaded new interval: " << t0 << " to " << t1 << " sec" << endl;
	dt = t1 - t0;

	if (has_h) {
		// bind raw cell arrays for the new timestep (no cell->point averaging)
		Hcell1 = dataset1->GetCellData()->GetArray(Hindex);
		Ucell1 = dataset1->GetCellData()->GetArray(Uindex);
		if (has_phi) Phicell1 = dataset1->GetCellData()->GetArray(Phiindex);
	}
	else if (cell2Pointdata) {
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


	if (recompute_betah_every_timstep) {
		// Compute height function betah
		if (input2d) {
			// Calculate betah for all points i LSgrid
			double z, welev, seabed, pNew[3];
			for (int k = 0; k < nl; k++) {
				for (int i = 0; i < nx; i++) {
					//cout << i << ", " << j << endl;
					dataset1->GetPoint(i, k, 0, pNew);
					z = pNew[1];
					dataset1->GetPoint(i, 0, 0, pNew);
					seabed = pNew[1];
					dataset1->GetPoint(i, nl - 1, 0, pNew);
					welev = pNew[1];
					betah[k * nx + i] = welev == seabed ? 0. : 1 + z2s(z, welev, -seabed);
					//cout << "seabed:" << seabed << " welev: " << welev << " betah:" << betah[k * nx  + i]<< endl;
				}
			}
		}
		else {
			//cout << "nx, ny, nl: " << nx << ", " << ny << ", " << nl << endl;
			// Calculate betah for all points i LSgrid
			double z, welev, seabed, pNew[3];
			for (int k = 0; k < nl; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						//cout << i << ", " << j << ", " << k << endl;
						dataset1->GetPoint(i, j, k, pNew);
						z = pNew[2];
						dataset1->GetPoint(i, j, 0, pNew);
						seabed = pNew[2];
						dataset1->GetPoint(i, j, nl - 1, pNew);
						welev = pNew[2];

						betah[k * ny * nx + j * nx + i] = welev == seabed ? 0. : 1. + z2s(z, welev, -seabed);
						//cout << "seabed:" << seabed << " welev: " << welev << " z: " << z << " betah:" << betah[k * ny * nx + j * nx + i] << endl;
					}
				}
			}
		}
	}
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
		selevs[k] = betah[k* ny * nx + nyp * nx + nxp];
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
double* VTKreader::trilinear_interpolation(double* res, double tpt, double xpt, double ypt, double zpt) {

	// array to store results, where 0 = eta, 1=seabed, 2=u, 3=v, 4=w

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

	// No point in calculating seabed twice, since seabed is fixed in time.
	/*dataset1->GetPoint(nxp, nyp, 0, pNew);
	D00 = pNew[2];
	dataset1->GetPoint(nxp, clip(nyp + 1, 0, ny - 1), 0, pNew);
	D01 = pNew[2];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), nyp, 0, pNew);
	D10 = pNew[2];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), clip(nyp + 1, 0, ny - 1), 0, pNew);
	D11 = pNew[2];*/

	C0 = C00 * (1. - xd) + C10 * xd;
	C1 = C01 * (1. - xd) + C11 * xd;
	//D0 = D00 * (1. - xd) + D10 * xd;
	//D1 = D01 * (1. - xd) + D11 * xd;

	double zb0 = C0 * (1. - yd) + C1 * yd;
	//double zb1 = D0 * (1. - yd) + D1 * yd;

	// set seabed elevation for point
	//res[1] =  zb0 * (1. - td) + zb1 * td;
	res[1] = zb0;

	//std::cout << "dy: " << dy << ", nyp: " << nyp << ", ypt: " << ypt << std::endl;

	double spt0 = wave_elev0 == zb0 ? -1. : std::max(z2s(std::min(zpt, wave_elev0), wave_elev0, -zb0), -1.);
	double spt1 = wave_elev1 == zb0 ? -1. : std::max(z2s(std::min(zpt, wave_elev1), wave_elev1, -zb0), -1.);
	

	int nsp0, nsp1;

	
	double sd0 = stretchInterpLocatorZ(spt0 + 1., &nsp0, nxp, nyp);
	double sd1 = stretchInterpLocatorZ(spt1 + 1., &nsp1, nxp, nyp);

	//cout << "sd0: " << sd0 << "z: " << zpt << endl;

	//int temp = clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1);
	
	//cout << "nxp: " << nxp << ", nyp: " << nyp << "nsp0: " << nsp0 << "sum: " << temp << endl;

	//exit(0);

	//cout << "nsp0: " << nsp0 << " " << clip(nsp0 + 1, 0, nl - 1) << endl;

	// Trilinear interpolation.
	double* temp;
	double C000[24];
	temp = U0->GetTuple3(nsp0 * ny * nx + nyp * nx + nxp);
	C000[0] = temp[0];
	C000[1] = temp[1];
	C000[2] = temp[2];
	temp = U0->GetTuple3(nsp0 * ny * nx + nyp * nx + clip(nxp + 1, 0, nx - 1));
	C000[3] = temp[0];
	C000[4] = temp[1];
	C000[5] = temp[2];
	temp = U0->GetTuple3(nsp0 * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + nxp);
	C000[6] = temp[0];
	C000[7] = temp[1];
	C000[8] = temp[2];
	temp = U0->GetTuple3(nsp0 * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + clip(nxp + 1, 0, nx - 1));
	C000[9] = temp[0];
	C000[10] = temp[1];
	C000[11] = temp[2];
	temp = U0->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * ny * nx + nyp * nx + nxp);
	C000[12] = temp[0];
	C000[13] = temp[1];
	C000[14] = temp[2];
	temp = U0->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * ny * nx + nyp * nx + clip(nxp + 1, 0, nx - 1));
	C000[15] = temp[0];
	C000[16] = temp[1];
	C000[17] = temp[2];
	temp = U0->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + nxp);
	C000[18] = temp[0];
	C000[19] = temp[1];
	C000[20] = temp[2];
	temp = U0->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + clip(nxp + 1, 0, nx - 1));
	C000[21] = temp[0];
	C000[22] = temp[1];
	C000[23] = temp[2];

	//cout << "ulow: " << C000[0] << " " << C000[3] << " " << C000[6] << " " << C000[9] << endl;
	//cout << "uhigh: " << C000[12] << " " << C000[15] << " " << C000[18] << " " << C000[21] << endl;

	// Trilinear interpolation.
	double D000[24];
	temp = U1->GetTuple3(nsp1 * ny * nx + nyp * nx + nxp);
	D000[0] = temp[0];
	D000[1] = temp[1];
	D000[2] = temp[2];
	temp = U1->GetTuple3(nsp1 * ny * nx + nyp * nx + clip(nxp + 1, 0, nx - 1));
	D000[3] = temp[0];
	D000[4] = temp[1];
	D000[5] = temp[2];
	temp = U1->GetTuple3(nsp1 * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + nxp);
	D000[6] = temp[0];
	D000[7] = temp[1];
	D000[8] = temp[2];
	temp = U1->GetTuple3(nsp1 * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + clip(nxp + 1, 0, nx - 1));
	D000[9] = temp[0];
	D000[10] = temp[1];
	D000[11] = temp[2];
	temp = U1->GetTuple3(clip(nsp1 + 1, 0, nl - 1) * ny * nx + nyp * nx + nxp);
	D000[12] = temp[0];
	D000[13] = temp[1];
	D000[14] = temp[2];
	temp = U1->GetTuple3(clip(nsp1 + 1, 0, nl - 1) * ny * nx + nyp * nx + clip(nxp + 1, 0, nx - 1));
	D000[15] = temp[0];
	D000[16] = temp[1];
	D000[17] = temp[2];
	temp = U1->GetTuple3(clip(nsp1 + 1, 0, nl - 1) * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + nxp);
	D000[18] = temp[0];
	D000[19] = temp[1];
	D000[20] = temp[2];
	temp = U1->GetTuple3(clip(nsp1 + 1, 0, nl - 1) * ny * nx + clip(nyp + 1, 0, ny - 1) * nx + clip(nxp + 1, 0, nx - 1));
	D000[21] = temp[0];
	D000[22] = temp[1];
	D000[23] = temp[2];

	//double sd0 = nsp0 - floor(nsp0);
	//double sd1 = nsp1 - floor(nsp1);

	double CC00[3], CC01[3], CC10[3], CC11[3];
	CC00[0] = C000[0] * (1. - xd) + C000[3] * xd;
	CC00[1] = C000[1] * (1. - xd) + C000[4] * xd;
	CC00[2] = C000[2] * (1. - xd) + C000[5] * xd;

	CC01[0] = C000[6] * (1. - xd) + C000[9] * xd;
	CC01[1] = C000[7] * (1. - xd) + C000[10] * xd;
	CC01[2] = C000[8] * (1. - xd) + C000[11] * xd;

	CC10[0] = C000[12] * (1. - xd) + C000[15] * xd;
	CC10[1] = C000[13] * (1. - xd) + C000[16] * xd;
	CC10[2] = C000[14] * (1. - xd) + C000[17] * xd;

	CC11[0] = C000[18] * (1. - xd) + C000[21] * xd;
	CC11[1] = C000[19] * (1. - xd) + C000[22] * xd;
	CC11[2] = C000[20] * (1. - xd) + C000[23] * xd;

	//cout << "ulow: " << CC00[0] << " " << CC01[0] << endl;
	//cout << "uhigh: " << CC10[0] << " " << CC11[0] << endl;

	//std::cout << int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0)) << ", " << ceil(nxp) << ", " << ceil(nyp) << ", " << ceil(nsp0) << std::endl;

	//std::cout << C010 << ", " << C110 << ", " << C011 << ", " << C111 << std::endl;

	double DD00[3], DD01[3], DD10[3], DD11[3];
	DD00[0] = D000[0] * (1. - xd) + D000[3] * xd;
	DD00[1] = D000[1] * (1. - xd) + D000[4] * xd;
	DD00[2] = D000[2] * (1. - xd) + D000[5] * xd;

	DD01[0] = D000[6] * (1. - xd) + D000[9] * xd;
	DD01[1] = D000[7] * (1. - xd) + D000[10] * xd;
	DD01[2] = D000[8] * (1. - xd) + D000[11] * xd;

	DD10[0] = D000[12] * (1. - xd) + D000[15] * xd;
	DD10[1] = D000[13] * (1. - xd) + D000[16] * xd;
	DD10[2] = D000[14] * (1. - xd) + D000[17] * xd;

	DD11[0] = D000[18] * (1. - xd) + D000[21] * xd;
	DD11[1] = D000[19] * (1. - xd) + D000[22] * xd;
	DD11[2] = D000[20] * (1. - xd) + D000[23] * xd;

	double CC0[3], CC1[3], DD0[3], DD1[3];
	CC0[0] = CC00[0] * (1. - yd) + CC01[0] * yd;
	CC0[1] = CC00[1] * (1. - yd) + CC01[1] * yd;
	CC0[2] = CC00[2] * (1. - yd) + CC01[2] * yd;
	CC1[0] = CC10[0] * (1. - yd) + CC11[0] * yd;
	CC1[1] = CC10[1] * (1. - yd) + CC11[1] * yd;
	CC1[2] = CC10[2] * (1. - yd) + CC11[2] * yd;

	//cout << "ulow: " << CC0[0]  << endl;
	//cout << "uhigh: " << CC1[0] << endl;

	DD0[0] = DD00[0] * (1. - yd) + DD01[0] * yd;
	DD0[1] = DD00[1] * (1. - yd) + DD01[1] * yd;
	DD0[2] = DD00[2] * (1. - yd) + DD01[2] * yd;
	DD1[0] = DD10[0] * (1. - yd) + DD11[0] * yd;
	DD1[1] = DD10[1] * (1. - yd) + DD11[1] * yd;
	DD1[2] = DD10[2] * (1. - yd) + DD11[2] * yd;

	
	// 2D data lives in the xy plane with y vertical; the velocity vector is
	// stored as (u_x, w, 0). Map to the 3D convention where z is the vertical
	// axis: VeloX = u_x (comp 0), VeloZ = w (comp 1), VeloY = 0 (no transverse
	// component exists in 2D data). This keeps the reader self-consistent with
	// the spatial convention above, where the query z is used as the vertical.
	res[2] = (CC0[0] * (1. - sd0) + CC1[0] * sd0) * (1. - td) + (DD0[0] * (1. - sd1) + DD1[0] * sd1) * td; // u (horizontal)
	res[4] = (CC0[1] * (1. - sd0) + CC1[1] * sd0) * (1. - td) + (DD0[1] * (1. - sd1) + DD1[1] * sd1) * td; // w (vertical, vts comp 1)
	res[3] = 0.0; // v (transverse — none in 2D data)

	//cout << "ulow:" << CC0[0] << ", uhigh: "<< CC1[0] << endl;

	return res;
}

/* Function for bilinear interpolation on a cartesian evenly spaced mesh in x direction and
 streched mesh in y direction*/
double* VTKreader::bilinear_interpolation_xy(double* res, double tpt, double xpt, double zpt) {

	// array to store results, where 0 = eta, 1=seabed, 2=u, 3=v, 4=w

	float nxp_temp;
	double xd = std::modf(clip((xpt - bounds[0]) / dx, 0., nx - 1.), &nxp_temp);

	int nxp = int(nxp_temp);

	double pNew[3];

	// get values of eta
	dataset0->GetPoint(nxp, nl - 1, 0, pNew);
	double C0 = pNew[1];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), nl - 1, 0, pNew);
	double C1 = pNew[1];

	dataset1->GetPoint(nxp, nl - 1, 0, pNew);
	double D0 = pNew[1];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), nl - 1, 0, pNew);
	double D1 = pNew[1];
	

	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double wave_elev0 = C0 * (1. - xd) + C1 * xd;
	double wave_elev1 = D0 * (1. - xd) + D1 * xd;

	// set seabed elevation for point
	res[0] = wave_elev0 * (1. - td) + wave_elev1 * td;

	// get values of seabed
	dataset0->GetPoint(nxp, 0, 0, pNew);
	C0 = pNew[1];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), 0, 0, pNew);
	C1 = pNew[1];
	

	dataset1->GetPoint(nxp, 0, 0, pNew);
	D0 = pNew[1];	
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), 0, 0, pNew);
	D1 = pNew[1];

	double zb0 = C0 * (1. - xd) + C1 * xd;
	double zb1 = D0 * (1. - xd) + D1 * xd;

	// set seabed elevation for point
	res[1] =  zb0 * (1. - td) + zb1 * td;


	//std::cout << "dy: " << dy << ", nyp: " << nyp << ", ypt: " << ypt << std::endl;

	double spt0 = wave_elev0 == zb0 ? -1. : std::max(z2s(std::min(zpt, wave_elev0), wave_elev0, -zb0), -1.);
	double spt1 = wave_elev1 == zb0 ? -1. : std::max(z2s(std::min(zpt, wave_elev1), wave_elev1, -zb1), -1.);
	

	int nsp0, nsp1;

	
	double sd0 = stretchInterpLocatorZ(spt0 + 1., &nsp0, nxp, 0);
	double sd1 = stretchInterpLocatorZ(spt1 + 1., &nsp1, nxp, 0);

	//int temp = clip(nxp + 1, 0, nx - 1) * ny * nl + clip(nyp + 1, 0, ny - 1) * nl + clip(nsp0 + 1, 0, nl - 1);
	//cout << "wave_elev: " << wave_elev0 << ", " << wave_elev1 << " sb: " << zb0 << ", " << zb1 << endl;
	//cout << "zpt: " << zpt <<" spt0: " << spt0 << endl;
	//cout << "nxp: " << nxp << ", " << "nsp0: " << nsp0 << endl;

	//exit(0);
	// Bilinear interpolation.
	double* temp;
	double CC00[8];
	temp = U0->GetTuple3(nsp0 * nx + nxp);
	//cout << "0u0:" << temp[0] << " u1:" << temp[1] << " u2:" << temp[2] << endl;
	CC00[0] = temp[0];
	CC00[4] = temp[1];
	temp = U0->GetTuple3(nsp0 * nx + clip(nxp + 1, 0, nx - 1));
	//cout << "1u0:" << temp[0] << " u1:" << temp[1] << " u2:" << temp[2] << endl;
	CC00[1] = temp[0];
	CC00[5] = temp[1];
	temp = U0->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * nx + nxp);
	//cout << "2u0:" << temp[0] << " u1:" << temp[1] << " u2:" << temp[2] << endl;
	CC00[2] = temp[0];
	CC00[6] = temp[1];
	temp = U0->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * nx + clip(nxp + 1, 0, nx - 1));
	//cout << "3u0:" << temp[0] << " u1:" << temp[1] << " u2:" << temp[2] << endl;
	CC00[3] = temp[0];
	CC00[7] = temp[1];

	//cout << "CC00:" << CC00[0] << " CC01:" << CC00[1] << " CC02:" << CC00[2] << " CC03:" << CC00[3] <<endl;
	//cout << "CC04:" << CC00[4] << " CC05:" << CC00[5] << " CC06:" << CC00[6] << " CC07:" << CC00[7] <<endl;

	

	double DD00[8];
	temp = U1->GetTuple3(nsp0 * nx + nxp);
	DD00[0] = temp[0];
	DD00[4] = temp[1];
	temp = U1->GetTuple3(nsp0 * nx + clip(nxp + 1, 0, nx - 1));
	DD00[1] = temp[0];
	DD00[5] = temp[1];
	temp = U1->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * nx + nxp);
	DD00[2] = temp[0];
	DD00[6] = temp[1];
	temp = U1->GetTuple3(clip(nsp0 + 1, 0, nl - 1) * nx + clip(nxp + 1, 0, nx - 1));
	DD00[3] = temp[0];
	DD00[7] = temp[1];

	//double sd0 = nsp0 - floor(nsp0);
	//double sd1 = nsp1 - floor(nsp1);

	//std::cout << int(ceil(nxp) * ny * nl + ceil(nyp) * nl + ceil(nsp0)) << ", " << ceil(nxp) << ", " << ceil(nyp) << ", " << ceil(nsp0) << std::endl;

	//std::cout << C010 << ", " << C110 << ", " << C011 << ", " << C111 << std::endl;

	double CC0[2], CC1[2], DD0[2], DD1[2];
	CC0[0] = CC00[0] * (1. - xd) + CC00[1] * xd;
	CC1[0] = CC00[2] * (1. - xd) + CC00[3] * xd;
	CC0[1] = CC00[4] * (1. - xd) + CC00[5] * xd;
	CC1[1] = CC00[6] * (1. - xd) + CC00[7] * xd;
	DD0[0] = DD00[0] * (1. - xd) + DD00[1] * xd;
	DD1[0] = DD00[2] * (1. - xd) + DD00[3] * xd;
	DD0[1] = DD00[4] * (1. - xd) + DD00[5] * xd;
	DD1[1] = DD00[6] * (1. - xd) + DD00[7] * xd;

	
	res[2] = (CC0[0] * (1. - sd0) + CC1[0] * sd0) * (1. - td) + (DD0[0] * (1. - sd1) + DD1[0] * sd1) * td; // u
	res[4] = (CC0[1] * (1. - sd0) + CC1[1] * sd0) * (1. - td) + (DD0[1] * (1. - sd1) + DD1[1] * sd1) * td; // w
	res[3] = 0.; // v

	return res;
}




/* Function for bi interpolation on a cartesian evenly spaced mesh*/
double* VTKreader::bilinear_interpolation(double* res, double tpt, double xpt, double ypt) {

	// array to store results, where 0 = eta, 1=seabed

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
	res[1] = zb0 * (1. - td) + zb1 * td;

	return res;
}


/* Function for bi interpolation on a cartesian evenly spaced mesh*/
double* VTKreader::linear_interpolation(double* res, double tpt, double xpt) {

	// array to store results, where 0 = eta, 1=seabed

	float nxp_temp;
	double xd = std::modf(clip((xpt - bounds[0]) / dx, 0., nx - 1.), &nxp_temp);
	int nxp = int(nxp_temp);
	
	double pNew[3];

	// get values of eta
	dataset0->GetPoint(nxp, nl - 1, 0, pNew);
	double C0 = pNew[1];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), nl - 1, 0, pNew);
	double C1 = pNew[1];
	

	dataset1->GetPoint(nxp, nl - 1, 0, pNew);
	double D0 = pNew[1];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), nl - 1, 0, pNew);
	double D1 = pNew[1];
	
	double td = std::min(1., std::max(0., (tpt - t0) / dt));

	double wave_elev0 = C0 * (1. - xd) + C1 * xd;
	double wave_elev1 = D0 * (1. - xd) + D1 * xd;

	// set seabed elevation for point
	res[0] = wave_elev0 * (1. - td) + wave_elev1 * td;

	// get values of seabed
	dataset0->GetPoint(nxp, 0, 0, pNew);
	C0 = pNew[1];
	dataset0->GetPoint(clip(nxp + 1, 0, nx - 1), 0, 0, pNew);
	C1 = pNew[1];

	dataset1->GetPoint(nxp, 0, 0, pNew);
	D0 = pNew[1];
	dataset1->GetPoint(clip(nxp + 1, 0, nx - 1), 0, 0, pNew);
	D1 = pNew[1];

	double zb0 = C0 * (1. - xd) + C1 * xd;
	double zb1 = D0 * (1. - xd) + D1 * xd;

	// set seabed elevation for point
	res[1] = zb0 * (1. - td) + zb1 * td;

	return res;
}

// ---------------------------------------------------------------------------
// Detect the optional cell-centred "h" (layer thickness) and "phi" (non-hydrostatic
// pressure) arrays in dataset1, and bind the raw cell-centred velocity/h/phi arrays
// for the multilayer reconstruction. The has_h/has_phi flags are sticky (checked once).
// ---------------------------------------------------------------------------
void VTKreader::locate_multilayer_fields() {
	vtkIdType nCell = dataset1->GetCellData()->GetNumberOfArrays();
	for (vtkIdType i = 0; i < nCell; i++) {
		const char* nm = dataset1->GetCellData()->GetArrayName(i);
		if (!nm) continue;
		if (strcmp(nm, "h") == 0)   { Hindex = i;   has_h = true; }
		if (strcmp(nm, "phi") == 0) { Phiindex = i; has_phi = true; }
	}
	if (has_h) {
		Hcell1 = dataset1->GetCellData()->GetArray(Hindex);
		// raw (un-averaged) cell velocity for the cell-centre reconstruction
		Ucell1 = dataset1->GetCellData()->GetArray(Uindex);
	}
	if (has_phi) Phicell1 = dataset1->GetCellData()->GetArray(Phiindex);
}

// ---------------------------------------------------------------------------
// h-based multilayer reconstruction. Fills res[0..5] = {eta, seabed, u, v, w, phi}.
// eta is the exact per-column sum of layer thicknesses h; velocity is reconstructed
// from the cell centres (PPR or linear); phi is reconstructed on the layer interfaces
// (phi_cell[k] lives on the lower interface of cell k; the surface value is 0).
// Horizontal blending follows the layers, i.e. at constant sigma (fraction of the
// local water column), consistent with the Lagrangian stretched grid.
// ---------------------------------------------------------------------------
double* VTKreader::multilayer_interpolation(double* res, double tpt, double xpt, double ypt, double zpt) {
	const int ncx = nx - 1;                       // cells in x
	const int ncy = input2d ? 1 : (ny - 1);       // cells in y
	const int ncell = nl - 1;                     // layers (cells in z)

	// locate cell column (interpolate between cell centres, which sit at (i+0.5)*d)
	auto locate = [](double q, double q0, double d, int nc, double& frac) -> int {
		double c = (q - q0) / d - 0.5;
		if (c < 0.0) { frac = 0.0; return 0; }
		if (c > nc - 1) { frac = 0.0; return std::max(0, nc - 1); }
		double fi; frac = std::modf(c, &fi);
		int idx = (int)fi;
		if (idx >= nc - 1) { idx = std::max(0, nc - 2); frac = (nc >= 2) ? 1.0 : 0.0; }
		return idx;
	};
	double xd = 0.0, yd = 0.0;
	int ic = locate(xpt, bounds[0], dx, ncx, xd);
	int jc = input2d ? 0 : locate(ypt, bounds[2], dy, ncy, yd);
	const double td = std::min(1.0, std::max(0.0, (tpt - t0) / dt));

	// vertical buffers, reused across calls (thread-local so it stays re-entrant)
	static thread_local std::vector<double> zi, zc, hh, uu, vv, ww, pI;
	if ((int)zi.size() != ncell + 1) {
		zi.resize(ncell + 1); pI.resize(ncell + 1);
		zc.resize(ncell); hh.resize(ncell); uu.resize(ncell); vv.resize(ncell); ww.resize(ncell);
	}

	// seabed (bathymetry) at a cell-column: average of the four (or two) bounding vertices
	auto seabed_col = [&](vtkStructuredGrid* ds, int icol, int jcol) -> double {
		icol = std::max(0, std::min(icol, ncx - 1));
		jcol = std::max(0, std::min(jcol, ncy - 1));
		double p[3];
		if (input2d) {
			ds->GetPoint(icol, 0, 0, p);                       double a = p[1];
			ds->GetPoint(std::min(icol + 1, nx - 1), 0, 0, p); double b = p[1];
			return 0.5 * (a + b);
		}
		double s = 0.0;
		for (int dj = 0; dj <= 1; dj++) for (int di = 0; di <= 1; di++) {
			ds->GetPoint(std::min(icol + di, nx - 1), std::min(jcol + dj, ny - 1), 0, p);
			s += p[2];
		}
		return 0.25 * s;
	};

	// lightweight build (pass 1): only eta = zb + sum(h) and zb, reading h alone.
	auto build_eta = [&](vtkStructuredGrid* ds, vtkDataArray* Ha, int icol, int jcol,
	                     double& eta_c, double& zb_c) {
		zb_c = seabed_col(ds, icol, jcol);
		icol = std::max(0, std::min(icol, ncx - 1));
		jcol = std::max(0, std::min(jcol, ncy - 1));
		double sumh = 0.0;
		for (int k = 0; k < ncell; k++) {
			long cid = input2d ? (long)k * ncx + icol
			                   : (long)k * ncy * ncx + (long)jcol * ncx + icol;
			sumh += Ha ? Ha->GetTuple1(cid) : 0.0;
		}
		eta_c = zb_c + sumh;
	};

	// build a full column (geometry + fields) for cell-column (icol,jcol) of dataset ds
	auto build = [&](vtkStructuredGrid* ds, vtkDataArray* Ua, vtkDataArray* Ha, vtkDataArray* Pa,
	                 int icol, int jcol, double& eta_c, double& zb_c) {
		zb_c = seabed_col(ds, icol, jcol);
		icol = std::max(0, std::min(icol, ncx - 1));
		jcol = std::max(0, std::min(jcol, ncy - 1));
		zi[0] = zb_c;
		for (int k = 0; k < ncell; k++) {
			long cid = input2d ? (long)k * ncx + icol
			                   : (long)k * ncy * ncx + (long)jcol * ncx + icol;
			hh[k] = Ha ? Ha->GetTuple1(cid) : 0.0;
			double* uvw = Ua->GetTuple3(cid);
			uu[k] = uvw[0]; vv[k] = uvw[1]; ww[k] = uvw[2];
			zi[k + 1] = zi[k] + hh[k];
			zc[k] = zi[k] + 0.5 * hh[k];
			pI[k] = Pa ? Pa->GetTuple1(cid) : 0.0;
		}
		pI[ncell] = 0.0;                          // phi = 0 at the free surface
		eta_c = zi[ncell];
	};

	// ---- pass 1: per-corner eta (sum of h) and seabed, then blended eta/zb ----
	double etaC[2][4], zbC[4];
	int ci[4] = { ic, std::min(ic + 1, ncx - 1), ic, std::min(ic + 1, ncx - 1) };
	int cj[4] = { jc, jc, std::min(jc + 1, ncy - 1), std::min(jc + 1, ncy - 1) };
	for (int c = 0; c < 4; c++) {
		double e, zb;
		build_eta(dataset0, Hcell0, ci[c], cj[c], e, zb); etaC[0][c] = e; zbC[c] = zb;
		build_eta(dataset1, Hcell1, ci[c], cj[c], e, zb); etaC[1][c] = e;
	}
	auto bilin = [&](double q00, double q10, double q01, double q11) {
		return (q00 * (1 - xd) + q10 * xd) * (1 - yd) + (q01 * (1 - xd) + q11 * xd) * yd;
	};
	double eta0 = bilin(etaC[0][0], etaC[0][1], etaC[0][2], etaC[0][3]);
	double eta1 = bilin(etaC[1][0], etaC[1][1], etaC[1][2], etaC[1][3]);
	double zbq  = bilin(zbC[0], zbC[1], zbC[2], zbC[3]);
	double eta_q = eta0 * (1 - td) + eta1 * td;
	res[0] = eta_q;
	res[1] = zbq;

	// sigma of the query in the (interpolated) local water column, clamped to [-1, 0]
	double s_q = (eta_q <= zbq) ? -1.0 : std::max(-1.0, std::min(0.0, z2s(std::min(zpt, eta_q), eta_q, -zbq)));

	// ---- pass 2: reconstruct u,v,w (cell centres) and phi (interfaces) per column ----
	double U[2][4], V[2][4], W[2][4], P[2][4];
	for (int dsn = 0; dsn < 2; dsn++) {
		vtkStructuredGrid* ds = dsn ? dataset1 : dataset0;
		vtkDataArray* Ua = dsn ? Ucell1 : Ucell0;
		vtkDataArray* Ha = dsn ? Hcell1 : Hcell0;
		vtkDataArray* Pa = dsn ? Phicell1 : Phicell0;
		for (int c = 0; c < 4; c++) {
			double e, zb;
			build(ds, Ua, Ha, Pa, ci[c], cj[c], e, zb);
			double z_local = e + s_q * (e - zb);      // same sigma -> this column's physical z
			U[dsn][c] = ppr_reconstruct(zi.data(), zc.data(), uu.data(), hh.data(), ncell, z_local, vertical_interp);
			// 2D data stores velocity as (u_x, w, 0): vertical is comp 1, no transverse.
			// 3D data stores (u_x, u_y, w): transverse is comp 1, vertical is comp 2.
			if (input2d) {
				V[dsn][c] = 0.0;
				W[dsn][c] = ppr_reconstruct(zi.data(), zc.data(), vv.data(), hh.data(), ncell, z_local, vertical_interp);
			} else {
				V[dsn][c] = ppr_reconstruct(zi.data(), zc.data(), vv.data(), hh.data(), ncell, z_local, vertical_interp);
				W[dsn][c] = ppr_reconstruct(zi.data(), zc.data(), ww.data(), hh.data(), ncell, z_local, vertical_interp);
			}
			// phi: piecewise-linear on the interface values (point data, not cell averages)
			if (has_phi) {
				double zz = std::max(zi[0], std::min(zi[ncell], z_local));
				int m = 0; while (m < ncell - 1 && zz > zi[m + 1]) m++;
				double wgt = (zi[m + 1] > zi[m]) ? (zz - zi[m]) / (zi[m + 1] - zi[m]) : 0.0;
				P[dsn][c] = pI[m] * (1 - wgt) + pI[m + 1] * wgt;
			} else {
				P[dsn][c] = 0.0;
			}
		}
	}
	double u0 = bilin(U[0][0], U[0][1], U[0][2], U[0][3]), u1 = bilin(U[1][0], U[1][1], U[1][2], U[1][3]);
	double v0 = bilin(V[0][0], V[0][1], V[0][2], V[0][3]), v1 = bilin(V[1][0], V[1][1], V[1][2], V[1][3]);
	double w0 = bilin(W[0][0], W[0][1], W[0][2], W[0][3]), w1 = bilin(W[1][0], W[1][1], W[1][2], W[1][3]);
	double p0 = bilin(P[0][0], P[0][1], P[0][2], P[0][3]), p1 = bilin(P[1][0], P[1][1], P[1][2], P[1][3]);
	res[2] = u0 * (1 - td) + u1 * td;
	res[3] = v0 * (1 - td) + v1 * td;
	res[4] = w0 * (1 - td) + w1 * td;
	res[5] = p0 * (1 - td) + p1 * td;
	return res;
}

double VTKreader::phi(double tpt, double xpt, double ypt, double zpt) {
	if (!has_phi) return 0.0;
	double r[6]; get_kinematics(tpt, xpt, ypt, zpt, r); return r[5];
}

bool VTKreader::CheckTime(double tpt) {
	/* Checks to see if the time tpt is within the interval t0 to t1. If so, returns true*/
	if ((tpt + t_start) > t1) {
		//std::cout << "t0: " << t0 << ", t1: " << t1 << ", tpt: " << tpt << std::endl;
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
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int m = 0; m < nl; m++) {
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
					double zpt1 = s2z(betah[i * ny * nl + j * nl + m] - 1., eta1_temp, -seabed);
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
					double zpt0 = s2z(betah[m * ny * nx + j * nx + i] - 1., eta0_temp, -seabed);
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

