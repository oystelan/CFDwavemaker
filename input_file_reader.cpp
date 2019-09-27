
#include "input_file_reader.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
//#include <cctype>
//#include <locale>
#include <algorithm>

using namespace std;


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
