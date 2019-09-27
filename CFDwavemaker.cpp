// This program has the soul purpose of providing wave kinematics input to any type
// of CFD program which can be linked up as a dynamic link library or statically.
// The program is created and updated by the Oeystein Lande. The link library
// compiles on both linux and windows. The following wave theory/types are currently
// supported:
// linear wave theory (2D and 3D), second order wave theory (2D and 3D) (sharma &
// dean), wave paddle theory (2D only at the moment)
//
//
// Current version: v2.0.1
// Date: 2019-09-30
// --------------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <ctime>
#include <cctype>
#include <locale>
#include "omp.h"
#include <algorithm>
#include <direct.h> // windows only functio
#include "CFDwavemaker.h" 

//#include <fftw3.h>

using namespace std;

double *UX;
double *UY;
double *UZ;
double *UXL;
double *UYL;
double *UZL;
double *ETA;

// Variables
int nfreq, ndir, wavetype, extmet, pertmet, meth, bandwidth, n_timesteps, rampswitch, normalizeA, spreadfunc;
int NX, NY, NZ, NXL, NYL, NZL;
double ampl, depth, s, mtheta, tofmax, fpoint[2], trampdata[3], xrampdata[3], yrampdata[3];
double fp, alpha_z, alpha_u, ramp_time, x_pos, y_pos, current_speed, wave_length, wave_height;



//fftw_plan p;

// Declaration of pointers where data will be stored
double* w;
double* Ampspec;
double* k;
double* thetaA;
double* D;
double* phas;
double* dsum2;
double* PD_time;
double* PD_ampl;
double* PD_velo;
double* PD_eta;
double* domainsize;
double* index;

extern "C" {
#include "Stokes5.h"
}

Stokes5 wave; // declare wave

int initialized;
int initsurf = 0;
int initkin = 0;
double dx, dy, dz, dxl, dyl, dzl;


//#define GetCurrentDir _getcwd

/*void testfft() {
	fftw_complex *in, *out;

	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 64);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 64);

	fftw_free(in); fftw_free(out);

}*/

void wait(int seconds)
{
	clock_t endwait;
	endwait = clock() + seconds * CLOCKS_PER_SEC;
	while (clock() < endwait) {}
}

// trim from start (in place)
static inline void ltrim(string& s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(), [](int ch) {
		return !isspace(ch);
		}));
}

// trim from end (in place)
static inline void rtrim(string& s) {
	s.erase(find_if(s.rbegin(), s.rend(), [](int ch) {
		return !isspace(ch);
		}).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(string& s) {
	ltrim(s);
	rtrim(s);
}

double sum(double ll[], int nsum) {
	double ss = 0.0;
	for (int i = 0; i < nsum; i++) {
		ss += ll[i];
	}
	return ss;
}

int check_license()
{

	int licensecheck = 1;
  // Date for program to stop working
  int expyear = 2020;
  int expmonth = 12;
  int expday = 31;

  time_t t = time(0);   // get time now

	#if defined(_MSC_VER)

      struct tm now;

      localtime_s(&now, &t);



    	cout << (now.tm_year + 1900) << '-'
    		<< (now.tm_mon + 1) << '-'
    		<< now.tm_mday
    		<< endl;
    	if ((now.tm_year + 1900) > expyear) {
    		licensecheck = 0;
    	}
    	else if ((now.tm_year + 1900) == expyear) {
    		if ((now.tm_mon + 1) > expmonth) {
    			licensecheck = 0;
    		}
    	}
    	else if ((now.tm_year + 1900) == expyear) {
    		if ((now.tm_mon + 1) == expmonth) {
    			if (now.tm_mday > expday) {
    				licensecheck = 0;
    			}
    		}
    	}
  #else
      struct tm * now;
      //time(&t)
      now = localtime(&t);

    	cout << (now->tm_year + 1900) << '-'
    		<< (now->tm_mon + 1) << '-'
    		<< now->tm_mday
    		<< endl;
    	if ((now->tm_year + 1900) > expyear) {
    		licensecheck = 0;
    	}
    	else if ((now->tm_year + 1900) == expyear) {
    		if ((now->tm_mon + 1) > expmonth) {
    			licensecheck = 0;
    		}
    	}
    	else if ((now->tm_year + 1900) == expyear) {
    		if ((now->tm_mon + 1) == expmonth) {
    			if (now->tm_mday > expday) {
    				licensecheck = 0;
    			}
    		}
    	}

  #endif
	return licensecheck;
}


//string GetCurrentWorkingDir(void) {
//	char buff[FILENAME_MAX];
//	GetCurrentDir(buff, FILENAME_MAX);
//	std::string current_working_dir(buff);
//	return current_working_dir;
//}




int read_inputdata_v2() {
	string lineA;
	ifstream fid;
	string res;

	// READ INPUT FILE AND REMOVE COMMENT LINES
	fid.open("./waveinput.dat");

	// check one step up in the folder tree (this is used in the latest comflow version)
	if (fid.fail()) {
		fid.open("../waveinput.dat");
	}
	// Error check
	if (fid.fail()) {
		cerr << "Could not open file (is it really there?) " << endl;
		return -1;
		exit(1);
	}
	else {
		cout << "Reading data from file: waveinput.dat..." << endl;
	}
	while (fid.good()) {
		getline(fid, lineA);
		//cout << lineA << endl;
		lineA.erase(find(lineA.begin(), lineA.end(), '#'), lineA.end());
		if (lineA.length() > 0) {
			res += lineA + "\n";

		}
	}
	fid.close();

	istringstream buf;
	istringstream f(res);
	//get and write data lines
	while (!f.eof()) {
		getline(f, lineA);
		trim(lineA);
		cout << lineA << endl;

		if (!lineA.compare("[wave type]")) {
			getline(f, lineA);
			wavetype = stoi(lineA);
		}
		if (!lineA.compare("[general wave data]")) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> ampl;
			buf >> normalizeA;
			buf >> depth;
			buf >> mtheta;
			buf.clear();
		}
		if (!lineA.compare("[perturbation]")) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> extmet;
			buf >> pertmet;
			buf >> bandwidth;
			buf.clear();
		}
		if (!lineA.compare("[wave origin]")) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> tofmax;
			buf >> fpoint[0];
			buf >> fpoint[1];
			buf.clear();
		}
		if (!lineA.compare("[ramps]")) {
			// read time ramp data
			getline(f, lineA);
			buf.str(lineA);
			buf >> ramp_time;
			buf.clear();
			if (ramp_time > 0) {
				rampswitch = 1;
			}
			else {
				rampswitch = 0;
			}
			// read ramp in x direction
			getline(f, lineA);
			buf.str(lineA);
			buf >> xrampdata[0];
			buf >> xrampdata[1];
			buf >> xrampdata[2];
			buf.clear();
			// read ramp in y direction
			getline(f, lineA);
			buf.str(lineA);
			buf >> yrampdata[0];
			buf >> yrampdata[1];
			buf >> yrampdata[2];
			buf.clear();
		}
		if (!lineA.compare("[wave type specific data]")) {
			// Wave specific data
			if (wavetype == 1) {
				// read spreading function data
				getline(f, lineA);
				buf.str(lineA);
				buf >> spreadfunc;
				buf >> s;
				buf.clear();

				// read number of frequencies and directions
				getline(f, lineA);
				buf.str(lineA);
				buf >> nfreq;
				buf >> ndir;
				buf.clear();

				// Read frequency data (omega, Sw and K)
				double* w_temp = new double[nfreq];
				double* Sw = new double[nfreq];
				double* k_temp = new double[nfreq];
				double* phas_temp = new double[nfreq];
				for (int i = 0; i < nfreq; i++) {
					getline(f, lineA);
					buf.str(lineA);
					buf >> w_temp[i];
					buf >> Sw[i];
					buf >> k_temp[i];
					buf >> phas_temp[i];
					buf.clear();
				}


				// Read directions
				double* theta_temp = new double[ndir];
				for (int i = 0; i < ndir; i++) {
					getline(f, lineA);
					buf.str(lineA);
					buf >> theta_temp[i];
					buf.clear();
				}

				// Normalize amplitude spectrum if Normalize is switched on
				double* Ampspec_temp = new double[nfreq];
				if (normalizeA) {
					for (int i = 0; i < nfreq; i++) {
						Ampspec_temp[i] = ampl * Sw[i] / sum(Sw, nfreq);
					}
				}
				else {
					for (int i = 0; i < nfreq; i++) {
						Ampspec_temp[i] = ampl * Sw[i];
					}
				}
				delete[] Sw;


				// Assign wave spreading based on specified spreading function parameters
				D = new double[nfreq * ndir];
				double dsum = 0.;
				dsum2 = new double[nfreq];
				if (spreadfunc == 0) { // Uniform-distribution
					for (int i = 0; i < ndir; i++) {
						D[i] = 1.0;
						dsum += D[i];
					}
					for (int i = 0; i < ndir; i++) {
						D[i] = D[i] / dsum;
					}
					// Repeat Directional distribution nfreq times
					for (int i = 0; i < nfreq; i++) {
						for (int j = 0; j < ndir; j++) {
							D[i * ndir + j] = D[j];
						}
					}
				}
				else if (spreadfunc == 1) { // cos(theta)^s
					for (int i = 0; i < ndir; i++) {
						D[i] = pow(cos((theta_temp[i] - (mtheta * PI / 180.))), s);
						dsum += D[i];
					}
					for (int i = 0; i < ndir; i++) {
						D[i] = D[i] / dsum;
					}
					// Repeat Directional distribution nfreq times
					for (int i = 0; i < nfreq; i++) {
						for (int j = 0; j < ndir; j++) {
							D[i * ndir + j] = D[j];
						}
					}
				}
				else if (spreadfunc == 2) { // cos(theta/2)^2s
					for (int i = 0; i < ndir; i++) {
						D[i] = pow(cos((theta_temp[i] - (mtheta * PI / 180.)) / 2.), 2.0 * s);
						dsum += D[i];
					}
					for (int i = 0; i < ndir; i++) {
						D[i] = D[i] / dsum;
					}
					// Repeat Directional distribution nfreq times
					for (int i = 0; i < nfreq; i++) {
						for (int j = 0; j < ndir; j++) {
							D[i * ndir + j] = D[j];
						}
					}
				}
				else if (spreadfunc == 3) { // Ewans simplified spreading function (frequency dependent spreading)
					for (int i = 0; i < nfreq; i++) {
						if (((w[i] / (2. * PI)) / fp) < 1.) {
							s = 15.5 * pow((w[i] / (2. * PI)) / fp, 9.47);
						}
						else {
							s = 13.1 * pow((w[i] / (2. * PI)) / fp, -1.94);
						}
						dsum2[i] = 0.;
						for (int j = 0; j < ndir; j++) {
							D[i * ndir + j] = pow(cos((theta_temp[j] - (mtheta * PI / 180.)) / 2.), 2. * s);
							dsum2[i] += D[i * ndir + j];
						}
					}
					for (int i = 0; i < nfreq; i++) {
						for (int j = 0; j < ndir; j++) {
							D[i * ndir + j] = D[i * ndir + j] / dsum2[i];
						}
					}
				}
				else if (spreadfunc == 4) { //johannessen cos(theta/2)^s special case (johannessen 1997)
					for (int i = 0; i < ndir; i++) {
						D[i] = pow(cos((theta_temp[i] - (mtheta * PI / 180.)) / 2.), s);
						dsum += D[i];
					}
					for (int i = 0; i < ndir; i++) {
						D[i] = D[i] / dsum;
					}
					// Repeat Directional distribution nfreq times
					for (int i = 0; i < nfreq; i++) {
						for (int j = 0; j < ndir; j++) {
							D[i * ndir + j] = D[j];
						}
					}
				}

				// Restack frequency and direction dimentions into 1 dimentional arrays
				w = new double[nfreq * ndir];
				k = new double[nfreq * ndir];
				phas = new double[nfreq * ndir];
				Ampspec = new double[nfreq * ndir];
				thetaA = new double[nfreq * ndir];

				for (int i = 0; i < nfreq; i++) {
					for (int j = 0; j < ndir; j++) {
						w[i * ndir + j] = w_temp[i];
						k[i * ndir + j] = k_temp[i];
						Ampspec[i * ndir + j] = Ampspec_temp[i];
						phas[i * ndir + j] = phas_temp[i];
						thetaA[i * ndir + j] = theta_temp[j];

					}
				}
				delete[] Ampspec_temp, w_temp, phas_temp, k_temp, theta_temp;

			}
		}
		if (!lineA.compare("[grid]")) {

		}
	}
	cout << "Input file read successfully." << endl;

}


int read_inputdata()
{
	string lineA;
	ifstream fid;
	string res;

	// READ INPUT FILE AND REMOVE COMMENT LINES
	fid.open("./waveinput.dat");

	// check one step up in the folder tree (this is used in the latest comflow version)
	if (fid.fail()) {
		fid.open("../waveinput.dat");
	}
	// Error check
	if (fid.fail()) {
		cerr << "Could not open file (is it really there?) " << endl;
		return -1;
		exit(1);
	}
	else {
		cout << "Reading data from file: waveinput.dat..." << endl;
	}
	while (fid.good()) {
		getline(fid, lineA);
		//cout << lineA << endl;
		lineA.erase(find(lineA.begin(), lineA.end(), '#'), lineA.end());
		if (lineA.length() > 0) {
			res += lineA + "\n";

		}
	}
	fid.close();
	//cout << res << endl;
	istringstream f(res);
	getline(f, lineA);

	wavetype = stoi(lineA);
	//cout << "The following wave type was chosen: " << wavetype << endl;

	istringstream buf;

	// ----------------------------------------------------------------------------------------------
	// WAVE TYPE 1 data read
	// The traditional method of specification (integrates both frequency and directions)
	// ----------------------------------------------------------------------------------------------
	if (wavetype == 1) {

		cout << endl << "WaveType: 1" << endl;
		cout << "Description: Irregular wave field with multiple directional components " << endl;
		cout << "-----------------------------------------------------------------------" << endl;

		// read Line 1
		//
		getline(f, lineA);
		buf.str(lineA);
		buf >> ampl;
		buf >> normalizeA;
		buf >> depth;
		buf >> mtheta;
		buf >> ramp_time;
		buf.clear();

		if (ramp_time > 0) {
			rampswitch = 1;
		}
		else {
			rampswitch = 0;
		}


		// read line 2
		getline(f, lineA);
		buf.str(lineA);
		buf >> extmet;
		buf >> pertmet;
		buf >> bandwidth;
		buf.clear();



		// read Line 3 % Phase adjustements: time of max, x-focuspoint, y-focuspoint
		getline(f, lineA);
		buf.str(lineA);
		buf >> tofmax;
		buf >> fpoint[0];
		buf >> fpoint[1];
		buf.clear();



		// read Line 5 % values for x ramp
		getline(f, lineA);
		buf.str(lineA);
		buf >> xrampdata[0];
		buf >> xrampdata[1];
		buf >> xrampdata[2];
		buf.clear();

		// read Line 6 % values for y ramp
		getline(f, lineA);
		buf.str(lineA);
		buf >> yrampdata[0];
		buf >> yrampdata[1];
		buf >> yrampdata[2];
		buf.clear();

		// read spreading function data
		getline(f, lineA);
		buf.str(lineA);
		buf >> spreadfunc;
		buf >> s;
		buf.clear();

		// read number of frequencies and directions
		getline(f, lineA);
		buf.str(lineA);
		buf >> nfreq;
		buf >> ndir;
		buf.clear();

		// Read frequency data (omega, Sw and K)
		double* w_temp = new double[nfreq];
		double* Sw = new double[nfreq];
		double* k_temp = new double[nfreq];
		double* phas_temp = new double[nfreq];
		for (int i = 0; i < nfreq; i++) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> w_temp[i];
			buf >> Sw[i];
			buf >> k_temp[i];
			buf >> phas_temp[i];
			buf.clear();
		}


		// Read directions
		double* theta_temp = new double[ndir];
		for (int i = 0; i < ndir; i++) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> theta_temp[i];
			buf.clear();
		}

		// Normalize amplitude spectrum if Normalize is switched on
		double* Ampspec_temp = new double[nfreq];
		if (normalizeA) {
			for (int i = 0; i < nfreq; i++) {
				Ampspec_temp[i] = ampl * Sw[i] / sum(Sw, nfreq);
			}
		}
		else {
			for (int i = 0; i < nfreq; i++) {
				Ampspec_temp[i] = ampl * Sw[i];
			}
		}
		delete[] Sw;


		// Assign wave spreading based on specified spreading function parameters
		D = new double[nfreq * ndir];
		double dsum = 0.;
		dsum2 = new double[nfreq];
		if (spreadfunc == 0) { // Uniform-distribution
			for (int i = 0; i < ndir; i++) {
				D[i] = 1.0;
				dsum += D[i];
			}
			for (int i = 0; i < ndir; i++) {
				D[i] = D[i] / dsum;
			}
			// Repeat Directional distribution nfreq times
			for (int i = 0; i < nfreq; i++) {
				for (int j = 0; j < ndir; j++) {
					D[i * ndir + j] = D[j];
				}
			}
		}
		else if (spreadfunc == 1) { // cos(theta)^s
			for (int i = 0; i < ndir; i++) {
				D[i] = pow(cos((theta_temp[i] - (mtheta * PI / 180.))), s);
				dsum += D[i];
			}
			for (int i = 0; i < ndir; i++) {
				D[i] = D[i] / dsum;
			}
			// Repeat Directional distribution nfreq times
			for (int i = 0; i < nfreq; i++) {
				for (int j = 0; j < ndir; j++) {
					D[i * ndir + j] = D[j];
				}
			}
		}
		else if (spreadfunc == 2) { // cos(theta/2)^2s
			for (int i = 0; i < ndir; i++) {
				D[i] = pow(cos((theta_temp[i] - (mtheta * PI / 180.)) / 2.), 2.0 * s);
				dsum += D[i];
			}
			for (int i = 0; i < ndir; i++) {
				D[i] = D[i] / dsum;
			}
			// Repeat Directional distribution nfreq times
			for (int i = 0; i < nfreq; i++) {
				for (int j = 0; j < ndir; j++) {
					D[i * ndir + j] = D[j];
				}
			}
		}
		else if (spreadfunc == 3) { // Ewans simplified spreading function (frequency dependent spreading)
			for (int i = 0; i < nfreq; i++) {
				if (((w[i] / (2. * PI)) / fp) < 1.) {
					s = 15.5 * pow((w[i] / (2. * PI)) / fp, 9.47);
				}
				else {
					s = 13.1 * pow((w[i] / (2. * PI)) / fp, -1.94);
				}
				dsum2[i] = 0.;
				for (int j = 0; j < ndir; j++) {
					D[i * ndir + j] = pow(cos((theta_temp[j] - (mtheta * PI / 180.)) / 2.), 2. * s);
					dsum2[i] += D[i * ndir + j];
				}
			}
			for (int i = 0; i < nfreq; i++) {
				for (int j = 0; j < ndir; j++) {
					D[i * ndir + j] = D[i * ndir + j] / dsum2[i];
				}
			}
		}
		else if (spreadfunc == 4) { //johannessen cos(theta/2)^s special case (johannessen 1997)
			for (int i = 0; i < ndir; i++) {
				D[i] = pow(cos((theta_temp[i] - (mtheta * PI / 180.)) / 2.), s);
				dsum += D[i];
			}
			for (int i = 0; i < ndir; i++) {
				D[i] = D[i] / dsum;
			}
			// Repeat Directional distribution nfreq times
			for (int i = 0; i < nfreq; i++) {
				for (int j = 0; j < ndir; j++) {
					D[i * ndir + j] = D[j];
				}
			}
		}

		// Restack frequency and direction dimentions into 1 dimentional arrays
		w = new double[nfreq * ndir];
		k = new double[nfreq * ndir];
		phas = new double[nfreq * ndir];
		Ampspec = new double[nfreq * ndir];
		thetaA = new double[nfreq * ndir];

		for (int i = 0; i < nfreq; i++) {
			for (int j = 0; j < ndir; j++) {
				w[i * ndir + j] = w_temp[i];
				k[i * ndir + j] = k_temp[i];
				Ampspec[i * ndir + j] = Ampspec_temp[i];
				phas[i * ndir + j] = phas_temp[i];
				thetaA[i * ndir + j] = theta_temp[j];

			}
		}
		delete[] Ampspec_temp, w_temp, phas_temp, k_temp, theta_temp;
		cout << "Read input wave file completed." << endl;
	}

	// ------------------------------------------------------------------------------------------------------------------
	// WAVE TYPE 2 data read (One directional component pr. frequency)
	// ------------------------------------------------------------------------------------------------------------------

	else if (wavetype == 2) {

		cout << endl << "WaveType: 2" << endl;
		cout << "Description: Irregular wave field with single directional component pr frequency " << endl;
		cout << "---------------------------------------------------------------------------------" << endl;
		// read Line 1
		getline(f, lineA);
		buf.str(lineA);
		buf >> ampl;
		buf >> normalizeA;
		buf >> depth;
		buf >> mtheta;
		buf >> ramp_time;
		buf.clear();
		if (ramp_time > 0) {
			rampswitch = 1;
		}
		else {
			rampswitch = 0;
		}

		// read Line 2
		getline(f, lineA);
		buf.str(lineA);
		buf >> extmet;
		buf >> pertmet;
		buf >> bandwidth;
		buf.clear();


		// read Line 3
		getline(f, lineA);
		buf.str(lineA);
		buf >> tofmax;
		buf >> fpoint[0];
		buf >> fpoint[1];
		buf.clear();

		// read Line 4
		getline(f, lineA);
		buf.str(lineA);
		buf >> xrampdata[0];
		buf >> xrampdata[1];
		buf >> xrampdata[2];
		buf.clear();

		// read Line 5
		getline(f, lineA);
		buf.str(lineA);
		buf >> yrampdata[0];
		buf >> yrampdata[1];
		buf >> yrampdata[2];
		buf.clear();

		// read number wave components
		getline(f, lineA);
		buf.str(lineA);
		buf >> nfreq;
		buf.clear();
		ndir = 1;



		// Read frequency data (omega, ampltude, wavenumber, phase (rad), direction (rad))
		w = new double[nfreq];
		double* Sw = new double[nfreq];
		k = new double[nfreq];
		phas = new double[nfreq];
		thetaA = new double[nfreq];
		D = new double[nfreq];
		for (int i = 0; i < nfreq; i++) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> w[i];
			buf >> Sw[i];
			buf >> k[i];
			buf >> phas[i];
			buf >> thetaA[i];
			buf.clear();
			D[i] = 1.0;
		}

		// Normalize amplitude spectrum if Normalize is switched on
		Ampspec = new double[nfreq];
		if (normalizeA) {
			for (int i = 0; i < nfreq; i++) {
				Ampspec[i] = ampl * Sw[i] / sum(Sw, nfreq);
			}
		}
		else {
			for (int i = 0; i < nfreq; i++) {
				Ampspec[i] = ampl * Sw[i];
			}
		}
		delete[] Sw;

		cout << "Wave type 2 read successfully." << endl;
	}

	// ------------------------------------------------------------------------------------------------------------------
	// WAVE TYPE 3 data read (One directional component pr. frequency)
	// ------------------------------------------------------------------------------------------------------------------
	// Description: same wave type as 2, but optimized for initializing a large basin, which requires some additional
	// parameters


	else if (wavetype == 3) {

		cout << endl << "WaveType: 3" << endl;
		cout << "Description: Irregular wave field with single directional component pr frequency " << endl;
		cout << "---------------------------------------------------------------------------------" << endl;
		// read Line 1
		getline(f, lineA);
		buf.str(lineA);
		buf >> ampl;
		buf >> normalizeA;
		buf >> depth;
		buf >> mtheta;
		buf >> ramp_time;
		buf.clear();
		if (ramp_time > 0) {
			rampswitch = 1;
		}
		else {
			rampswitch = 0;
		}

		// read Line 2
		getline(f, lineA);
		buf.str(lineA);
		buf >> extmet;
		buf >> pertmet;
		buf >> bandwidth;
		buf.clear();


		// read Line 3
		getline(f, lineA);
		buf.str(lineA);
		buf >> tofmax;
		buf >> fpoint[0];
		buf >> fpoint[1];
		buf.clear();

		// read Line 4 Spatial ramp in x direction
		getline(f, lineA);
		buf.str(lineA);
		buf >> xrampdata[0];
		buf >> xrampdata[1];
		buf >> xrampdata[2];
		buf.clear();

		// read Line 5 Spatial ramp in y direction
		getline(f, lineA);
		buf.str(lineA);
		buf >> yrampdata[0];
		buf >> yrampdata[1];
		buf >> yrampdata[2];
		buf.clear();

		// read Line 6 Domain size
		domainsize = new double[7];
		getline(f, lineA);
		buf.str(lineA);
		buf >> domainsize[0];
		buf >> domainsize[1];
		buf >> domainsize[2];
		buf >> domainsize[3];
		buf >> domainsize[4];
		buf >> domainsize[5];
		buf >> domainsize[6];
		buf.clear();

		// read Line 7
		getline(f, lineA);
		buf.str(lineA);
		buf >> NX;
		buf >> NY;
		buf >> NZ;
		buf.clear();

		// read Line 8
		getline(f, lineA);
		buf.str(lineA);
		buf >> NXL;
		buf >> NYL;
		buf >> NZL;
		buf.clear();

		// read number wave components
		getline(f, lineA);
		buf.str(lineA);
		buf >> nfreq;
		buf.clear();
		ndir = 1;



		// Read frequency data (omega, ampltude, wavenumber, phase (rad), direction (rad))
		w = new double[nfreq];
		double* Sw = new double[nfreq];
		k = new double[nfreq];
		phas = new double[nfreq];
		thetaA = new double[nfreq];
		D = new double[nfreq];
		for (int i = 0; i < nfreq; i++) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> w[i];
			buf >> Sw[i];
			buf >> k[i];
			buf >> phas[i];
			buf >> thetaA[i];
			buf.clear();
			D[i] = 1.0;
		}

		// Normalize amplitude spectrum if Normalize is switched on
		Ampspec = new double[nfreq];
		if (normalizeA) {
			for (int i = 0; i < nfreq; i++) {
				Ampspec[i] = ampl * Sw[i] / sum(Sw, nfreq);
			}
		}
		else {
			for (int i = 0; i < nfreq; i++) {
				Ampspec[i] = ampl * Sw[i];
			}
		}
		delete[] Sw;

		cout << "Wave type 3 read successfully." << endl;
	}


	// ---------------------------------------------------------------------
	// PISTON TYPE WAVE MAKER
	// ---------------------------------------------------------------------
	else if (wavetype == 4) {
		cout << endl << "WaveType: 4" << endl;
		cout << "Description: Piston Type wave maker theory. Uses piston flap signal to calculate kinematics" << endl;
		cout << "---------------------------------------" << endl;
		// read alpha values
		getline(f, lineA);
		buf.str(lineA);
		buf >> alpha_z;
		buf >> alpha_u;
		buf.clear();

		getline(f, lineA);
		buf.str(lineA);
		buf >> n_timesteps;
		//n_timesteps = stoi(lineA);
		cout << "Number of timesteps: " << wavetype << endl;

		// declare some vectors to store piston data
		PD_time = new double[n_timesteps];
		PD_ampl = new double[n_timesteps];
		PD_velo = new double[n_timesteps];
		PD_eta = new double[n_timesteps];

		for (int i = 0; i < n_timesteps; i++) {
			getline(f, lineA);
			buf.str(lineA);
			buf >> PD_time[i];
			buf >> PD_ampl[i];
			buf >> PD_velo[i];
			buf >> PD_eta[i];
			buf.clear();
		}

	}

	// ----------------------------------------------------------------
	// WAVE TYPE 5 - Stokes 5th order
	// ----------------------------------------------------------------
	else if (wavetype == 5) {
		cout << endl << "WaveType: 5" << endl;
		cout << "Description: Stokes 5th order wave " << endl;
		cout << "---------------------------------------" << endl;
		// read Line 1
		getline(f, lineA);
		buf.str(lineA);
		buf >> wave_length;
		buf >> wave_height;
		buf >> depth;
		buf >> current_speed;
		buf >> mtheta;
		buf >> ramp_time;
		buf.clear();
		if (ramp_time > 0) {
			rampswitch = 1;
		}
		else {
			rampswitch = 0;
		}

		// read Line 2
		getline(f, lineA);
		buf.str(lineA);
		buf >> x_pos; // initial position of max
		buf >> y_pos; // position of still water level
		buf.clear();

		// read Line 4
		getline(f, lineA);
		buf.str(lineA);
		buf >> xrampdata[0];
		buf >> xrampdata[1];
		buf >> xrampdata[2];
		buf.clear();

		// read Line 5
		getline(f, lineA);
		buf.str(lineA);
		buf >> yrampdata[0];
		buf >> yrampdata[1];
		buf >> yrampdata[2];
		buf.clear();


		// set the properties of the wave
		set_stokes5_properties(&wave, wave_length, wave_height, current_speed, depth, G, x_pos, y_pos);
	}

	// Define which combination of methods to be used

	if (wavetype < 3) {
		// linear waves, expenential above swl
		if (extmet == 0 && pertmet == 0) {
			meth = 1;
			cout << "Linear Theory used. Exponential profile continued above the free surface" << endl;
		}
		// linear waves, constant value above swl
		else if (extmet == 1 && pertmet == 0) {
			meth = 2;
			cout << "Linear Theory used. Constant value (z=0) continued above the free surface" << endl;
		}
		// Second order solutions - Exponential
		else if (extmet == 0 && pertmet == 2) {
			meth = 3;
			cout << "Second order Theory used. Exponential profile continued above the free surface (unaccurate). Recommend using extmet==2" << endl;
		}

		// Second order solutions - Constant above Z
		else if (extmet == 1 && pertmet == 2) {
			meth = 4;
			cout << "Second order Theory used. Constant value (z=0) continued above the free surface. Recommend using extmet==2" << endl;

		}
		// Second order solutions + Taylor expansion
		else if (extmet == 2 && pertmet == 2) {
			meth = 5;
			cout << "Second order Theory used. Velocities extrapolated using Taylor expansion above free surface" << endl;

		}

	}
	// Special case where surface elevation and velocities are stored to a coarse grid and waves are
	else if (wavetype == 3) {
		meth = 8;
	}
	// Wave piston special case
	else if (wavetype == 4) {
		meth = 6;
	}
	// Stokes 5th wave
	else if (wavetype == 5) {
		meth = 7;
	}
	else {
		cerr << "What the? Illigal combination of methods!!! " << endl;
	}
	return 0;
}



//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
static double ramp(double x,double xsign, double xstart, double xend) {

	if (xsign > 0.){
		if (x <= xstart) {
			return 1.0;
		}
		else if (x >= xend) {
			return 0.0;
		}
		else {
			return 1. - ((x - xstart) / (xend - xstart));
		}
	}
	else {
		return 1.0;
	}
}


//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
static double timeramp(double t, double tsign, double tstart, double tend) {

	if (tsign > 0.) {
		if (t <= tstart) {
			return 0.0;
		}
		else if (t >= tend) {
			return 1.0;
		}
		else {
			return ((t - tstart) / (tend - tstart));
		}
	}
	else {
		return 1.0;
	}
}


// Linear interpolation function

int findNearestNeighbourIndex(double value, double *x, int len)
{
	double dist;
	int idx;
	int i;
	idx = -1;
	dist = DBL_MAX;
	for (i = 0; i < len; i++) {
		double newDist = value - x[i];
		if (newDist >= 0 && newDist < dist) {
			dist = newDist;
			idx = i;
		}
	}
	return idx;
}
double interp1(double *x, int x_tam, double *y, double xx)
{
	double dx, dy, slope, intercept, yy;
	int indiceEnVector;


	indiceEnVector = findNearestNeighbourIndex(xx, x, x_tam);
	if (indiceEnVector != -1) {
		dx = x[indiceEnVector + 1] - x[indiceEnVector];
		dy = y[indiceEnVector + 1] - y[indiceEnVector];
		slope = dy / dx;
		intercept = y[indiceEnVector] - x[indiceEnVector] * slope;
		yy = slope * xx + intercept;
	}
	else
		yy = DBL_MAX;

	return yy;
}


//--------------------------------------------------------------------
// -------------------------------------------------------------------
//
//  WAVEFUNCTIONS
//
//



/* wave elevation for a sinus wave */
double waveelev(double t, double xx, double yy) {

	double welev = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {

		phi = w[i] * tofmax+phas[i];
		welev += Ampspec[i] * D[i] *cos(k[i] * (cos(thetaA[i]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);

		//welev += A*cos(w[i] * t - k[i] * (cos(thetaA[i])*(xx - fpoint[0]) + sin(thetaA[i])*(yy - fpoint[1])) + phi);
	}


	return welev;

}


double waveelev_2order(double t, double xx, double yy) {



	//double eta1_t = 0;
	double eta2_t = 0;

	for (int i = 0; i < nfreq-1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));
				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
				}
				double alpha_nm_minus = (((w[ci] / w[cm]) + (w[cm] / w[ci])) + (G*G / (w[ci] * w[cm]))*
					((D_nm_minus - k[ci] * k[cm] * (gamma_nm + tanh(k[ci] * depth)*tanh(k[cm] * depth))) / (w[ci] * w[cm])));
				double alpha_nm_plus = (((w[ci] / w[cm]) + (w[cm] / w[ci])) + (G*G / (w[ci] * w[cm]))*
					((D_nm_plus - k[ci] * k[cm] * (gamma_nm - tanh(k[ci] * depth)*tanh(k[cm] * depth))) / (w[ci] * w[cm])));

				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				eta2_t += ((Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm])/(2.*G)) * (alpha_nm_minus*cos(phi_i - phi_m) + alpha_nm_plus*cos(phi_i + phi_m));



			}
		}
	}

	return eta2_t;

}


/* Second order horizontal velocity U for a sinus wave */
double uu_2order(double t, double xx, double yy, double zz) {



	//double eta1_t = 0;
	double usum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				usum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*cos(phi_i - phi_m)*(k[ci]*cos(thetaA[ci]+ (mtheta*PI / 180.))-k[cm]*cos(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_minus*(zz+depth))/cosh(k_nm_minus*depth)) +
					beta_nm_plus*cos(phi_i + phi_m)*(k[ci] * cos(thetaA[ci]+ (mtheta*PI / 180.)) + k[cm] * cos(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



			}
		}
	}
	return usum2;
}


/* Second order horizontal velocity V for a sinus wave */
double vv_2order(double t, double xx, double yy, double zz) {



	//double eta1_t = 0;
	double vsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}
				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				vsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*cos(phi_i - phi_m)*(k[ci] * sin(thetaA[ci]+ (mtheta*PI / 180.)) - k[cm] * sin(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_minus*(zz + depth)) / cosh(k_nm_minus*depth)) +
					beta_nm_plus*cos(phi_i + phi_m)*(k[ci] * sin(thetaA[ci]+ (mtheta*PI / 180.)) + k[cm] * sin(thetaA[cm]+ (mtheta*PI / 180.)))*(cosh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



			}
		}
	}
	return vsum2;
}



/* Second order vertical velocity component W for a sinus wave */
double ww_2order(double t, double xx, double yy, double zz) {

	//double eta1_t = 0;
	double wsum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adjusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));

				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm]+ (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm]+ (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				wsum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*k_nm_minus*sin(phi_i - phi_m)*(sinh(k_nm_minus*(zz + depth)) / cosh(k_nm_minus*depth)) +
					beta_nm_plus*k_nm_plus*sin(phi_i + phi_m)*(sinh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



			}
		}
	}
	return wsum2;
}

/* Second order component of velocity potential 
Warning: Unsure if this is finished...should be checked thoroughly at some point*/
double phi_pot_2order(double t, double xx, double yy, double zz) {
	//double eta1_t = 0;
	double phisum2 = 0;

	for (int i = 0; i < nfreq - 1; i++) {
		for (int j = 0; j < ndir; j++) {
			// Second order
			int ci = i*ndir + j;
			double phi_i = k[ci] * (cos(thetaA[ci] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[ci] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[ci] * t + w[ci] * tofmax + phas[ci];
			double Rn = k[ci] * tanh(k[ci] * depth);
			//cout << "Rn: " << Rn << endl;

			//// Adiusting Bandwidth for 2 order cut-off
			//if (i + 1 + f_bw<nfreqs) {
			//	p = i + 1 + f_bw;
			//}
			//else { p = nfreqs; }
			for (int m = i + 1; m < min(nfreq, i + bandwidth); m++) {
				int cm = m*ndir + j;
				double gamma_nm = cos(thetaA[ci] - thetaA[cm]);
				double k_nm_plus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] + (2.*k[ci] * k[cm] * gamma_nm));
				double k_nm_minus = sqrt(k[ci] * k[ci] + k[cm] * k[cm] - (2.*k[ci] * k[cm] * gamma_nm));

				double Rm = k[cm] * tanh(k[cm] * depth);


				double D_nm_plus = (sqrt(Rn) + sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) + sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) + sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm - Rn*Rm) /
					(pow(sqrt(Rn) + sqrt(Rm), 2.) - k_nm_plus*tanh(k_nm_plus*depth));
				double D_nm_minus = (sqrt(Rn) - sqrt(Rm))*(sqrt(Rm)*(k[ci] * k[ci] - Rn*Rn) - sqrt(Rn)*(k[cm] * k[cm] - Rm*Rm)) + 2.*pow(sqrt(Rn) - sqrt(Rm), 2.)*(k[ci] * k[cm] * gamma_nm + Rn*Rm) /
					(pow(sqrt(Rn) - sqrt(Rm), 2.) - k_nm_minus*tanh(k_nm_minus*depth));


				double beta_nm_minus = D_nm_minus / (2 * k[ci] * k[cm] * (w[ci] - w[cm]));
				double beta_nm_plus = D_nm_plus / (2 * k[ci] * k[cm] * (w[ci] + w[cm]));

				// Catch NaN when two equal frequency components interact
				if (w[ci] == w[cm]) {
					D_nm_minus = 0.;
					beta_nm_minus = 0.;
				}

				double phi_m = k[cm] * (cos(thetaA[cm] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[cm] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[cm] * t + w[cm] * tofmax + phas[cm];

				phisum2 += (Ampspec[ci] * D[ci] * Ampspec[cm] * D[cm] * w[ci] * w[cm]) *
					(beta_nm_minus*cos(phi_i - phi_m)*(k[ci] * cos(thetaA[ci] + (mtheta*PI / 180.)) - k[cm] * cos(thetaA[cm] + (mtheta*PI / 180.)))*(cosh(k_nm_minus*(zz + depth)) / cosh(k_nm_minus*depth)) +
						beta_nm_plus*cos(phi_i + phi_m)*(k[ci] * cos(thetaA[ci] + (mtheta*PI / 180.)) + k[cm] * cos(thetaA[cm] + (mtheta*PI / 180.)))*(cosh(k_nm_plus*(zz + depth)) / cosh(k_nm_plus*depth)));



			}
		}
	}
	return phisum2;
}

// First order velocity potential
double phi_pot(double t, double xx, double yy, double zz) {

	double phisum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		phisum += Ampspec[i] * D[i] * (w[i] / k[i]) * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return phisum;

}

double uu(double t, double xx, double yy, double zz) {

	double usum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		usum += cos(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return usum;

}


/* horizontal velocity U for a sinus wave */
double vv(double t, double xx, double yy, double zz) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		vsum += sin(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * (cosh(k[i] * (zz + depth)) / sinh(k[i] * depth))*cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return vsum;

}


/* vertical velocity for a sinus wave */
double ww(double t, double xx, double yy, double zz) {

	double wsum = 0.0;
	double phi;



	for (int i = 0; i< ndir*nfreq; i++) {

		phi = w[i] * tofmax + phas[i];
		wsum += D[i] * Ampspec[i] * w[i] * (sinh(k[i] * (zz + depth)) / sinh(k[i] * depth))*sin(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return wsum;

}

double pp(double t, double xx, double yy, double zz) {

	double psum = 0.0;
	double phi;


	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		psum += Ampspec[i] * D[i] * RHO * G * (cosh(k[i] * (zz + depth)) / cosh(k[i] * depth))*cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}

	return psum;

}

/* vertical velocity gradient at z=0 for velocity component U */
double phi_dxdz(double t, double xx, double yy) {

	double usum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		usum += cos(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * k[i] * cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}
	return usum;

}



/* vertical velocity gradient at z=0 for velocity component V */
double phi_dydz(double t, double xx, double yy) {

	double vsum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {
		phi = w[i] * tofmax + phas[i];
		vsum += sin(thetaA[i] + (mtheta*PI / 180.))* Ampspec[i] * D[i] * w[i] * k[i] *cos(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}
	return vsum;

}


/* vertical velocity gradient at z=0 for velocity component W */
double phi_dzdz(double t, double xx, double yy) {

	double wsum = 0.0;
	double phi;

	for (int i = 0; i< ndir*nfreq; i++) {

		phi = w[i] * tofmax + phas[i];
		wsum += D[i] * Ampspec[i] * w[i] * k[i] * (cosh(k[i] * depth) / sinh(k[i] * depth))*sin(k[i] * (cos(thetaA[i] + (mtheta*PI / 180.))*(xx - fpoint[0]) + sin(thetaA[i] + (mtheta*PI / 180.))*(yy - fpoint[1])) - w[i] * t + phi);
	}
	return wsum;

}


/* Horizontal velocity taken directly from the timeseries*/
double u_piston(double t) {
	double ux = interp1(PD_time, n_timesteps, PD_velo, t);

	return ux+alpha_u*ux;
}

/* Wave elevation taken directly from piston timeseries*/
double wave_elev_piston(double t) {
	return alpha_z*interp1(PD_time, n_timesteps, PD_eta, t);
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void initialize_kinematics(double tpt) {
	// Allocating memory for storage of surface elevation and velocities
	UX = new double[NX*NY*NZ];
	UY = new double[NX*NY*NZ];
	UZ = new double[NX*NY*NZ];

	UXL = new double[NXL*NYL*NZL];
	UYL = new double[NXL*NYL*NZL];
	UZL = new double[NXL*NYL*NZL];

	cout << "Memory allocation successful for storage of kinematics." << endl;

	dx = (domainsize[1] - domainsize[0]) / double(NX-1);
	dy = (domainsize[3] - domainsize[2]) / double(NY-1);
	dz = (domainsize[6] - domainsize[5]) / double(NZ-1);
	
	double dd = omp_get_wtime();

	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel // start parallell initialization
	{
		#pragma omp master
		cout << "Number of available threads: " << omp_get_num_threads() << endl;

		double xpt, ypt, zpt;
		double eta_temp;

		// Main grid
		#pragma omp for
		for (int i = 0; i < NX; i++) {
			xpt = domainsize[0] + dx*i;
			for (int j = 0; j < NY; j++) {
				ypt = domainsize[2] + dy*j;
				eta_temp = waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt);

				double Ux0 = uu(tpt, xpt, ypt, 0.0) + uu_2order(tpt, xpt, ypt, 0.0);
				double Uy0 = vv(tpt, xpt, ypt, 0.0) + vv_2order(tpt, xpt, ypt, 0.0);
				double Uz0 = ww(tpt, xpt, ypt, 0.0) + ww_2order(tpt, xpt, ypt, 0.0);

				double PHI_dxdz = phi_dxdz(tpt, xpt, ypt);
				double PHI_dydz = phi_dydz(tpt, xpt, ypt);
				double PHI_dzdz = phi_dzdz(tpt, xpt, ypt);

				for (int m = 0; m < NZ; m++) {
					zpt = domainsize[5] + dz*m;
					if (zpt > (eta_temp + dz)) {
						UX[i*NY*NZ + j*NZ + m] = 0.0;
						UY[i*NY*NZ + j*NZ + m] = 0.0;
						UZ[i*NY*NZ + j*NZ + m] = 0.0;
					}
					else if (zpt > 0.) {
						UX[i*NY*NZ + j*NZ + m] = Ux0 + PHI_dxdz*zpt;
						UY[i*NY*NZ + j*NZ + m] = Uy0 + PHI_dydz*zpt;
						UZ[i*NY*NZ + j*NZ + m] = Uz0 + PHI_dzdz*zpt;
					}
					else {
						UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt) + uu_2order(tpt, xpt, ypt, zpt);
						UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt) + vv_2order(tpt, xpt, ypt, zpt);
						UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt);
					}
					/*UX[i*NY*NZ + j*NZ + m] = uu(tpt, xpt, ypt, zpt);
					UY[i*NY*NZ + j*NZ + m] = vv(tpt, xpt, ypt, zpt);
					UZ[i*NY*NZ + j*NZ + m] = ww(tpt, xpt, ypt, zpt);*/
				}
			}
		}
	} // End parallel initialization

	cout << "Generation of upper domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	cout << "Initialization time: " << dd << " seconds." << endl;

	dd = omp_get_wtime();
	#pragma omp parallel // start parallel initialization
	{
		double xpt, ypt, zpt;
		// Secondary grid (coarse res at depth)
		dxl = (domainsize[1] - domainsize[0]) / double(NXL - 1);
		dyl = (domainsize[3] - domainsize[2]) / double(NYL - 1);
		dzl = (domainsize[5] - domainsize[4]) / double(NZL - 1);
		#pragma omp for
		for (int i = 0; i < NXL; i++) {
			xpt = domainsize[0] + dxl*i;
			for (int j = 0; j < NYL; j++) {
				ypt = domainsize[2] + dyl*j;
				for (int m = 0; m < NZL; m++) {
					zpt = domainsize[4] + dzl*m;
					UXL[i*NYL*NZL + j*NZL + m] = uu(tpt, xpt, ypt, zpt) + uu_2order(tpt, xpt, ypt, zpt);
					UYL[i*NYL*NZL + j*NZL + m] = vv(tpt, xpt, ypt, zpt) + vv_2order(tpt, xpt, ypt, zpt);
					UZL[i*NYL*NZL + j*NZL + m] = ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt);
				}
			}
		}
	} // End parallel
	cout << "Generation of lower domain kinematics data completed. ";
	dd = omp_get_wtime() - dd;
	cout << "Initialization time: " << dd << " seconds." << endl;
	cout << "Interpolation can commence..." << endl;
	initkin = 1;
}


// Precalculate velocityfield and surface elevation on coarse grid in case of WAVE TYPE 3
void initialize_surface_elevation(double tpt) {

	// Allocating memory for storage of surface elevation and velocities
	ETA = new double[NX*NY];

	cout << "Memory allocation successful for Surface elevation storage." << endl;

	dx = (domainsize[1] - domainsize[0]) / double(NX - 1);
	dy = (domainsize[3] - domainsize[2]) / double(NY - 1);
	
	double dd = omp_get_wtime();
	//omp_set_num_threads(1);
	omp_set_num_threads(omp_get_max_threads());

	#pragma omp parallel
	{
		double xpt, ypt;
		// Main grid
		#pragma omp for
		for (int i = 0; i < NX; i++) {
			xpt = domainsize[0] + dx*i;
			for (int j = 0; j < NY; j++) {
				ypt = domainsize[2] + dy*j;
				ETA[i*NY + j] = waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt);
			}
		}
	}
	dd = omp_get_wtime()-dd;

	cout << "Surface Elevation generated successfully. ";
	cout << "Initialization time: " << dd << " seconds." << endl;
	initsurf = 1;
}

/* Function for trilinear interpolation on a cartesian evenly spaced mesh*/
double trilinear_interpolation(double *VAR, double xpt, double ypt, double zpt) {
	double nxp = min(double(NX), max(0., (xpt - domainsize[0]) / dx));
	double nyp = min(double(NY), max(0., (ypt - domainsize[2]) / dy));
	double nzp = min(double(NZ), max(0., (zpt - domainsize[5]) / dz));

	double C000 = VAR[int(floor(nxp)*NY*NZ + floor(nyp)*NZ + floor(nzp))];
	double C001 = VAR[int(floor(nxp)*NY*NZ + floor(nyp)*NZ + ceil(nzp))];
	double C010 = VAR[int(floor(nxp)*NY*NZ + ceil(nyp)*NZ + floor(nzp))];
	double C011 = VAR[int(floor(nxp)*NY*NZ + ceil(nyp)*NZ + ceil(nzp))];
	double C100 = VAR[int(ceil(nxp)*NY*NZ + floor(nyp)*NZ + floor(nzp))];
	double C101 = VAR[int(ceil(nxp)*NY*NZ + floor(nyp)*NZ + ceil(nzp))];
	double C110 = VAR[int(ceil(nxp)*NY*NZ + ceil(nyp)*NZ + floor(nzp))];
	double C111 = VAR[int(ceil(nxp)*NY*NZ + ceil(nyp)*NZ + ceil(nzp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double zd = nzp - floor(nzp);

	double C00 = C000*(1. - xd) + C100*xd;
	double C01 = C001*(1. - xd) + C101*xd;
	double C10 = C010*(1. - xd) + C110*xd;
	double C11 = C011*(1. - xd) + C111*xd;

	double C0 = C00*(1. - yd) + C10*yd;
	double C1 = C01*(1. - yd) + C11*yd;

	return C0*(1. - zd) + C1*zd;
}
/* Function for trilinear interpolation on a cartesian evenly spaced mesh on the lower part of the domain*/
double trilinear_interpolationL(double *VAR, double xpt, double ypt, double zpt) {
	
	double nxp = min(double(NXL),max(0.,(xpt - domainsize[0]) / dxl));
	double nyp = min(double(NYL), max(0., (ypt - domainsize[2]) / dyl));
	double nzp = min(double(NZL), max(0., (zpt - domainsize[4]) / dzl));

	double C000 = VAR[int(floor(nxp)*NYL*NZL + floor(nyp)*NZL + floor(nzp))];
	double C001 = VAR[int(floor(nxp)*NYL*NZL + floor(nyp)*NZL + ceil(nzp))];
	double C010 = VAR[int(floor(nxp)*NYL*NZL + ceil(nyp)*NZL + floor(nzp))];
	double C011 = VAR[int(floor(nxp)*NYL*NZL + ceil(nyp)*NZL + ceil(nzp))];
	double C100 = VAR[int(ceil(nxp)*NYL*NZL + floor(nyp)*NZL + floor(nzp))];
	double C101 = VAR[int(ceil(nxp)*NYL*NZL + floor(nyp)*NZL + ceil(nzp))];
	double C110 = VAR[int(ceil(nxp)*NYL*NZL + ceil(nyp)*NZL + floor(nzp))];
	double C111 = VAR[int(ceil(nxp)*NYL*NZL + ceil(nyp)*NZL + ceil(nzp))];
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);
	double zd = nzp - floor(nzp);

	double C00 = C000*(1. - xd) + C100*xd;
	double C01 = C001*(1. - xd) + C101*xd;
	double C10 = C010*(1. - xd) + C110*xd;
	double C11 = C011*(1. - xd) + C111*xd;

	double C0 = C00*(1. - yd) + C10*yd;
	double C1 = C01*(1. - yd) + C11*yd;

	return C0*(1. - zd) + C1*zd;
}

/* bilinear interpolation function used to interpolate surface values on a regular evenly spaced grid*/
double bilinear_interpolation(double *VAR, double xpt, double ypt) {
	
	double nxp = min(double(NX), max(0.,(xpt - domainsize[0]) / dx));
	double nyp = min(double(NY), max(0.,(ypt - domainsize[2]) / dy));

	double C00 = VAR[int(floor(nxp)*NY + floor(nyp))];
	double C01 = VAR[int(floor(nxp)*NY + ceil(nyp))];
	double C10 = VAR[int(ceil(nxp)*NY + floor(nyp))];
	double C11 = VAR[int(ceil(nxp)*NY + ceil(nyp))];
	
	double xd = nxp - floor(nxp);
	double yd = nyp - floor(nyp);

	double C0 = C00*(1. - xd) + C10*xd;
	double C1 = C01*(1. - xd) + C11*xd;

	return C0*(1. - yd) + C1*yd;
}



//extern "C" {

// For comflow 3.9 pointer are sent. For later releases, the entire variables are sent

//EXPORT double VelocityX(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
double wave_VeloX(double xpt, double ypt, double zpt, double tpt)
{

	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
	// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*uu(tpt, xpt, ypt, zpt)*timeramp(tpt,rampswitch,0.,ramp_time);
	// Linear wave theory, constant profile used above free surface
	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*uu(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
	// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (uu(tpt, xpt, ypt, zpt) + uu_2order(tpt, xpt, ypt, zpt) )*timeramp(tpt, rampswitch, 0., ramp_time);
	// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * ( uu(tpt, xpt, ypt, z) + uu_2order(tpt, xpt, ypt, z) )*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (uu(tpt, xpt, ypt, z) + uu_2order(tpt, xpt, ypt, z) )*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((uu(tpt, xpt, ypt, z) + uu_2order(tpt, xpt, ypt, z)) + phi_dxdz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return u_piston(tpt);
	case 7:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*Stokes5_u(&wave, tpt, xpt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
	case 8:
		if (initkin == 0) {
			cout << "Generating kinematics for interpolation:" << endl;
			initialize_kinematics(0.0);
		}
		/*
		// Check if coordinates are within bounding box
		if (xpt < domainsize[0] || xpt > domainsize[1]) {
			cerr << "xpt: " << xpt << " out of bounds! Please extend interpolation box boundaries in x-direction" << endl;
		}
		else if (ypt < domainsize[2] || ypt > domainsize[3]) {
			cerr << "ypt: " << ypt << " out of bounds! Please extend interpolation box boundaries in y-direction" << endl;
		}
		else if (zpt < domainsize[4] || zpt > domainsize[6]) {
			cerr << "zpt: " << zpt << " out of bounds! Please extend interpolation box boundaries in z-direction" << endl;
		}
		*/
		if (zpt < domainsize[5]) {
			return trilinear_interpolationL(UXL, xpt, ypt, zpt);
		}
		else {
			return trilinear_interpolation(UX, xpt, ypt, zpt);
		}
	default:
		return 0.0;
	}


}

//
//EXPORT double VelocityY(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
double wave_VeloY(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*vv(tpt, xpt, ypt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface
	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*vv(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * ( vv(tpt, xpt, ypt, zpt) + vv_2order(tpt, xpt, ypt, zpt) )*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * ( vv(tpt, xpt, ypt, z) + vv_2order(tpt, xpt, ypt, z) )*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (vv(tpt, xpt, ypt, z) + vv_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((vv(tpt, xpt, ypt, z) + vv_2order(tpt, xpt, ypt, z)) + phi_dydz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return 0.0;
	case 7:
		return 0.0; // todo: implement directionality of stokes5th wave
	case 8:
		if (zpt < domainsize[5]) {
			return trilinear_interpolationL(UYL, xpt, ypt, zpt);
		}
		else {
			return trilinear_interpolation(UY, xpt, ypt, zpt);
		}
	default:
		return 0.0;
	}
}



//
//EXPORT double VelocityZ(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
double wave_VeloZ(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*ww(tpt, xpt, ypt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface
	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*ww(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt))*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*(ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z)) + phi_dzdz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return 0.0;
	case 7:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*Stokes5_v(&wave, tpt, xpt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
	case 8:
		if (zpt < domainsize[5]) {
			return trilinear_interpolationL(UZL, xpt, ypt, zpt);
		}
		else {
			return trilinear_interpolation(UZ, xpt, ypt, zpt);
		}
	default:
		return 0.0;
	}
}

//
double wave_DynPres(double xpt, double ypt, double zpt, double tpt)
//EXPORT double DynamicPressure(int& ii, int& jj, int& kk, double& xpt, double& ypt, double& zpt, double& tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = max(-depth, zpt);
	double z;

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*pp(tpt, xpt, ypt, zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface

	case 2:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*pp(tpt, xpt, ypt, z)*timeramp(tpt, rampswitch, 0., ramp_time);
	/*
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, zpt) + ww_2order(tpt, xpt, ypt, zpt))*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		z = min(0., zpt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*(ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		z = min(0., zpt);
		if (zpt <= 0) {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z))*timeramp(tpt, rampswitch, 0., ramp_time);
		}
		else {
			return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) *
				((ww(tpt, xpt, ypt, z) + ww_2order(tpt, xpt, ypt, z)) + phi_dzdz(tpt, xpt, ypt)*zpt)*timeramp(tpt, rampswitch, 0., ramp_time);
		}
	case 6:
		return 0.0;
		*/
	default:
		//cout << "WARNING: HIGHER ORDER PRESSURE CALCULATION NOT YET SUPPORTED." << endl;
		return 0.0;
	}
}

//
//EXPORT double SurfaceElevation(int& ii, int& jj, double& xpt, double& ypt, double& tpt)
double wave_SurfElev(double xpt, double ypt, double tpt)
{

	switch (meth) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		//return waveelev(tpt, xpt, ypt);
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*waveelev(tpt, xpt, ypt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Linear wave theory, constant profile used above free surface
	case 2:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*waveelev(tpt, xpt, ypt)*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, exponential profile used above free surface
	case 3:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt))*timeramp(tpt, rampswitch, 0., ramp_time);
		// Second order wave theory, constant profile used above free surface
	case 4:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 5:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2])) * (waveelev(tpt, xpt, ypt) + waveelev_2order(tpt, xpt, ypt))*timeramp(tpt, rampswitch, 0., ramp_time);
	case 6:
		return wave_elev_piston(tpt);
	case 7:
		return min(ramp(xpt, xrampdata[0], xrampdata[1], xrampdata[2]), ramp(ypt, yrampdata[0], yrampdata[1], yrampdata[2]))*Stokes5_eta(&wave, tpt, xpt)*timeramp(tpt, rampswitch, 0., ramp_time);
	case 8:
		if (initsurf == 0) {
			cout << "Initializing surface elevation storage:" << endl;
			initialize_surface_elevation(0.0);
			return bilinear_interpolation(ETA, xpt, ypt);
		}
		else {
			//cout << "asking for surface elevation..." << endl;
			//cout << xpt << " " << ypt << " " << tpt << endl;
			return bilinear_interpolation(ETA, xpt, ypt);
		}
	default:
		return 0.0;
	}
}



double wave_VFrac(double xpt, double ypt, double zpt, double tpt, double delta_cell)

	// Function which calculates volume fraction
	// assumes still water level z = 0
{
	double wwelev = wave_SurfElev(xpt,ypt,tpt);

	//cout << wwelev << endl;
	if (wwelev < zpt - (delta_cell / 2.)) {
		return 0.0;
	}
	else if (wwelev > zpt + (delta_cell / 2.)) {
		return 1.0;
	}
	else {
		// Calculate volume fraction for the given cell with size delta_cell and position zpt
		return (wwelev - (zpt - (delta_cell / 2.))) / delta_cell;
	}
	//return 0;
}




//EXPORT int Init(double& tmin_in, double& tmax_in)
int wave_Initialize()
{
	cout << "---------------------------------------" << endl;
	cout << "CFD WAVEMAKER v.1.07" << endl;
	cout << "---------------------------------------" << endl;
	
	// Check if license has expired
	if (check_license() == 1){
		cout << "License checks out...carry on..." << endl;
	}
	else {
		cout << "License for CFDwavemaker has expired. Please contact Oeystein Lande to update the program." << endl << endl;
		cout << "This program will auto-distruct in \n5..." << endl;
		wait(1);
		cout << "4..." << endl;
		wait(1);
		cout << "3..." << endl;
		wait(1);
		cout << "2..." << endl;
		wait(1);
		cout << "1..." << endl;
		wait(1);
		cout << "Bang!" << endl;

		return -1;
	}
	//for (int i = 0; i < nfreq; i++) {
	//	k[i] = pow(2. * pi * f[i], 2.) / 9.81;
	//}
	int i = read_inputdata();

	return 0;
}


int wave_Cleanup()
{
	if (wavetype == 1) {
		delete[] w, Ampspec, D, k, thetaA, phas;
	}
	else if (wavetype == 2) {
		delete[] w, Ampspec, k, thetaA, phas;
	}
	else if (wavetype == 3) {
		delete[] w, Ampspec, k, thetaA, phas;
		delete[] ETA, UX, UY, UZ, index;
	}
	else if (wavetype == 4) {
		delete[] PD_time, PD_ampl, PD_velo, PD_eta;
	}
	return 0;
}
//}



// external functions used by COMFLOW
double VelocityX(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	//cout << xpt << " " << ypt << " " << zpt << " " << wave_VeloX(xpt, ypt, zpt, time) << endl;
	return wave_VeloX(xpt, ypt, zpt, time);
}
double VelocityY(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return wave_VeloY(xpt, ypt, zpt, time);
}
double VelocityZ(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return wave_VeloZ(xpt, ypt, zpt, time);
}
double DynamicPressure(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return wave_DynPres(xpt, ypt, zpt, time);
}
double SurfaceElevation(int i, int j, double xpt, double ypt, double time) {
	return wave_SurfElev(xpt, ypt, time);
}
double VolumeFraction(double xpt, double ypt, double zpt, double time, double delta_xyz) {
	return wave_VFrac(xpt, ypt, zpt, time, delta_xyz);
}
double LiquidFillRatio(int i, int j, int kk, double xpt, double ypt, double zpt, double time) {
	return 0.0;
}
int Init(void* fptr) {
	return wave_Initialize();
}
int Prepare(double time, int mode) {
	return 0;
}
int Cleanup() {
	return wave_Cleanup();
}


int main() {
	//cout << GetCurrentWorkingDir() << endl;
	read_inputdata_v2();
}
