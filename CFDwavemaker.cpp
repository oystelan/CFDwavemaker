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
//#include <cstdlib>
//#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <ctime>
#include <algorithm>
//#include <direct.h> // windows only function
#//include <cctype>
//#include <locale>
#include "CFDwavemaker.h" 
#include "Stokes5.h"
#include "Irregular.h"
#include "Utils.h"
#include "Wavemaker.h"
#include "sgrid.h"

//#include <fftw3.h>


// Variables
//int nfreq, ndir, wavetype, extmet, pertmet, meth, bandwidth, n_timesteps, rampswitch, normalizeA, spreadfunc;


//double ampl, depth, s, mtheta, tofmax, fpoint[2], trampdata[3], xrampdata[3], yrampdata[3];
double x_pos, y_pos, tofmax, current_speed, wave_length, wave_height;
double mtheta;
double swl = 0.;

// Stokes 5 class
Stokes5 stokes5;

// Irregular class
Irregular irregular;

// Wavemaker theory class
Wavemaker wavemaker;

// Grid class
Grid gridclass;

// Grid class
sGrid sgrid;

// Ramp class
Ramp ramp;

//fftw_plan p;

int wavetype;

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
static inline void ltrim(std::string& s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(), [](int ch) {
		return !isspace(ch);
		}));
}

// trim from end (in place)
static inline void rtrim(std::string& s) {
	s.erase(find_if(s.rbegin(), s.rend(), [](int ch) {
		return !isspace(ch);
		}).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string& s) {
	ltrim(s);
	rtrim(s);
}

int numparams(std::string str)
{
	// breaking input into word using string stream 
	std::stringstream s(str); // Used for breaking words 
	std::string word; // to store individual words 

	int count = 0;
	while (s >> word)
		count++;
	return count;
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



    	std::cout << (now.tm_year + 1900) << '-'
    		<< (now.tm_mon + 1) << '-'
    		<< now.tm_mday
    		<< std::endl;
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

      std::cout << (now->tm_year + 1900) << '-'
    		<< (now->tm_mon + 1) << '-'
    		<< now->tm_mday
    		<< std::endl;
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
	std::string lineA;
	std::ifstream fid;
	std::string res;

	// READ INPUT FILE AND REMOVE COMMENT LINES
	fid.open("./waveinput.dat");

	// check one step up in the folder tree (this is used in the latest comflow version)
	if (fid.fail()) {
		fid.open("../waveinput.dat");
	}
	// Error check
	if (fid.fail()) {
		std::cerr << "Could not open file (is it really there?) " << std::endl;
		return -1;
		exit(1);
	}
	else {
		std::cout << "Reading data from file: waveinput.dat..." << std::endl;
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

	std::istringstream buf;
	std::istringstream f(res);
	//get and write data lines
	while (!f.eof()) {
		getline(f, lineA);
		trim(lineA);
		std::cout << lineA << std::endl;
		// Convension for numbering:
		// irregular wave theory variants: 1-10
		// wavemaker theory variants: 11-20
		// Regular wave theories: 21-30
		// HOSM and other: 31-40
		if (!lineA.compare("[wave type]")) {
			getline(f, lineA);
			trim(lineA);
			std::cout << lineA << std::endl;
			// check if valid wave type is given
			if (!lineA.compare("irregular") || !lineA.compare("1")) {
				wavetype = 1;
				std::cout << "Irregular perturbation wave theory specified" << std::endl;
			}
			else if (!lineA.compare("pistonwavemaker")) {
				wavetype = 11;
				std::cout << "Piston wave maker theory specified" << std::endl;
			}
			else if (!lineA.compare("spectral wave")) {
				wavetype = 31;
				std::cout << "spectral wave (HOSM) specified" << std::endl;
			}
			else if (!lineA.compare("stokes5")) {
				wavetype = 21;
				std::cout << "Regular 5th order Stokes wave specified" << std::endl;
			}
			else {
				std::cout << "Unknown wave type specified. Valid alternatives are:" << std::endl;
				std::cout << "irregular" << std::endl;
				std::cout << "pistonwavemaker" << std::endl;
				std::cout << "spectral wave" << std::endl;
				std::cout << "stokes5" << std::endl;
				//exit(1);
			}
		}
		if (!lineA.compare("[general input data]")) { //mandatory
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> depth;
			buf >> mtheta;
			buf.clear();
		}
		if (!lineA.compare("[normalize]")) { //optional
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> irregular.ampl;
			buf >> irregular.normalize;
			buf.clear();
			
		}
		if (!lineA.compare("[perturbation method]")) { //optional
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> irregular.extrapolation_met;
			buf >> irregular.order;
			buf >> irregular.dw_bandwidth;
			buf.clear();
		}
		if (!lineA.compare("[wave reference point]")) { //optional
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> tofmax;
			buf >> x_pos;
			buf >> y_pos;
			buf.clear();
		}
		if (!lineA.compare("[ramps]")) { //optional
			ramp.ramp_init = true;
			// read time ramp data
			
			// read time rampup
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> ramp.ramp_init_time_up;
			buf >> ramp.time_rampup_start;
			buf >> ramp.time_rampup_end;;
			buf.clear();
			// read time rampdown
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> ramp.ramp_init_time_down;
			buf >> ramp.time_rampdown_start;
			buf >> ramp.time_rampdown_end;;
			buf.clear();
			// read x-direction rampup
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> ramp.ramp_init_x_up;
			buf >> ramp.x_rampup_start;
			buf >> ramp.x_rampup_end;;
			buf.clear();
			// read x-direction rampdown
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> ramp.ramp_init_x_down;
			buf >> ramp.x_rampdown_start;
			buf >> ramp.x_rampdown_end;;
			buf.clear();
			// read y-direction rampup
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> ramp.ramp_init_y_up;
			buf >> ramp.y_rampup_start;
			buf >> ramp.y_rampup_end;;
			buf.clear();
			// read y-direction rampdown
			getline(f, lineA);
			std::cout << lineA << std::endl;
			buf.str(lineA);
			buf >> ramp.ramp_init_y_down;
			buf >> ramp.y_rampdown_start;
			buf >> ramp.y_rampdown_end;;
			buf.clear();
			
		}
		// Wave properties: this is where the wave type specific data is given
		if (!lineA.compare("[wave properties]")) {
			if (wavetype == 1) {
				// In case of irregular wave is specified
				getline(f, lineA);
				trim(lineA);
				std::cout << lineA << std::endl;
				// User defined wave. List of frequency components given
				// frequency, spectral ampl, wave number, phase, direction (rad)
				if (!lineA.compare("userdefined1")) {
					std::cout << "Irregular seas, one directional component for each frequency specified" << std::endl;
					// read number wave components
					getline(f, lineA);
					std::cout << "Number of frequency components: " << lineA << std::endl;
					buf.str(lineA);
					buf >> irregular.nfreq;
					buf.clear();
					irregular.ndir = 1;

					std::cout << "# OMEGA[rad / s]    A[m]           K             Phase[rad]     theta[rad]" << std::endl;
					double temp;
					for (int i = 0; i < irregular.nfreq; i++) {
						getline(f, lineA);
						std::cout << lineA << std::endl;
						buf.str(lineA);
						buf >> temp;
						irregular.omega.push_back(temp);
						buf >> temp;
						irregular.A.push_back(temp);
						buf >> temp;
						irregular.k.push_back(temp);
						buf >> temp;
						irregular.phase.push_back(temp);
						buf >> temp;
						irregular.theta.push_back(temp);
						buf.clear();
					}
				}
				// The traditional way of specifing frequency and direction as separate components S(f,theta) = S(f)*D(theta)
				else if (!lineA.compare("userdefined")) {
					std::cout << "Irregular seas, directional spreading defined separately" << std::endl;
					// read number of frequencies and directions
					getline(f, lineA);
					buf.str(lineA);
					buf >> irregular.nfreq;
					buf >> irregular.ndir;
					buf.clear();
					std::cout << "Number of frequency components: " << irregular.nfreq << std::endl;
					std::cout << "Number of directional components: " << irregular.ndir << std::endl;

					// Read frequency data (omega, Sw and K)
					double* omega_temp = new double[irregular.nfreq];
					double* Ampspec_temp = new double[irregular.nfreq];
					double* k_temp = new double[irregular.nfreq];
					double* phas_temp = new double[irregular.nfreq];
					std::cout << "# OMEGA[rad / s]    A[m]           K             Phase[rad]" << std::endl;
					for (int i = 0; i < irregular.nfreq; i++) {
						getline(f, lineA);
						std::cout << lineA << std::endl;
						buf.str(lineA);
						buf >> omega_temp[i];
						buf >> Ampspec_temp[i];
						buf >> k_temp[i];
						buf >> phas_temp[i];
						buf.clear();
					}

					// Read directional data
					double* theta_temp = new double[irregular.ndir];
					double* D_ampl_temp = new double[irregular.ndir];
					std::cout << "# Theta [rad] D(theta)" << std::endl;
					for (int i = 0; i < irregular.ndir; i++) {
						getline(f, lineA);
						buf.str(lineA);
						buf >> theta_temp[i];
						buf >> D_ampl_temp[i];
						buf.clear();
					}


					for (int i = 0; i < irregular.nfreq; i++) {
						for (int j = 0; j < irregular.ndir; j++) {
							irregular.omega.push_back(omega_temp[i]);
							irregular.k.push_back(k_temp[i]);
							irregular.A.push_back(Ampspec_temp[i] * D_ampl_temp[j]);
							irregular.phase.push_back(phas_temp[i]);
							irregular.theta.push_back(theta_temp[j]);

						}
					}
					delete[] Ampspec_temp, omega_temp, phas_temp, k_temp, theta_temp, D_ampl_temp;
				}

				
			}
			else if (wavetype == 11) {
				// Wavemaker theory (piston)
				// read alpha values
				getline(f, lineA);
				buf.str(lineA);
				buf >> wavemaker.alpha_z;
				buf >> wavemaker.alpha_u;
				buf.clear();

				getline(f, lineA);
				buf.str(lineA);
				buf >> wavemaker.n_timesteps;
				//n_timesteps = stoi(lineA);
				std::cout << "Number of timesteps: " << wavetype << std::endl;

				// declare some vectors to store piston data
				wavemaker.PD_time = new double[wavemaker.n_timesteps];
				wavemaker.PD_ampl = new double[wavemaker.n_timesteps];
				wavemaker.PD_velo = new double[wavemaker.n_timesteps];
				wavemaker.PD_eta = new double[wavemaker.n_timesteps];

				for (int i = 0; i < wavemaker.n_timesteps; i++) {
					getline(f, lineA);
					buf.str(lineA);
					buf >> wavemaker.PD_time[i];
					buf >> wavemaker.PD_ampl[i];
					buf >> wavemaker.PD_velo[i];
					buf >> wavemaker.PD_eta[i];
					buf.clear();
				}
			}
			else if (wavetype == 21) {
				// read Line 2
				getline(f, lineA);
				buf.str(lineA);
				buf >> stokes5.wave_length; // wave length
				buf >> stokes5.wave_height;  // Wave height
				buf.clear();

			}
			else {
				std::cout << "Unknown wavetype specified" << std::endl;
			}
		}

		if (!lineA.compare("[current speed]")) { //optional
			if (wavetype == 21){ // Stokes 5th
				getline(f, lineA);
				std::cout << lineA << std::endl;
				std::cout << "Current speed in m/s. Assumed inline with wave propagation direction." << std::endl;
				buf.str(lineA);
				buf >> stokes5.current;
				buf.clear();
			}
			else {
				std::cout << "Current not supported for selected wave type. Parameter ignored." << std::endl;
			}
		}
		if (!lineA.compare("[still water level]")) { //optional
			getline(f, lineA);
			std::cout << "still water level: " << lineA << "m" << std::endl;
			
			buf.str(lineA);
			buf >> swl;
			buf.clear();
		}

		if (!lineA.compare("[grid interpolation]")) {
			// Special case for fast initialization of domain
			getline(f, lineA);
			buf.str(lineA);
			buf >> gridclass.numgrids;
			buf.clear();

			for (int i = 0; i < gridclass.numgrids; i++) {
				getline(f, lineA);
				trim(lineA);
				std::cout << "Grid interpolation type: " << lineA << std::endl;
				if (!lineA.compare("init domain interp")) {
					getline(f, lineA);
					buf.str(lineA);
					buf >> gridclass.domain_start[0];
					buf >> gridclass.domain_end[0];
					buf >> gridclass.domain_start[1];
					buf >> gridclass.domain_end[1];
					buf >> gridclass.domain_start_L[2];
					buf >> gridclass.domain_end_L[2];
					buf >> gridclass.domain_end[2];
					gridclass.domain_start[2] = gridclass.domain_end_L[2];
					gridclass.domain_start_L[0] = gridclass.domain_start[0];
					gridclass.domain_start_L[1] = gridclass.domain_start[1];
					gridclass.domain_end_L[0] = gridclass.domain_end[0];
					gridclass.domain_end_L[1] = gridclass.domain_end[1];

					buf.clear();
					getline(f, lineA);
					buf.str(lineA);
					buf >> gridclass.NX;
					buf >> gridclass.NY;
					buf >> gridclass.NZ;
					buf.clear();
					getline(f, lineA);
					buf.str(lineA);
					buf >> gridclass.NXL;
					buf >> gridclass.NYL;
					buf >> gridclass.NZL;
					buf.clear();
					if (wavetype == 1)
						wavetype = 2;
				}
				if (!lineA.compare("wallx")) {
					getline(f, lineA);
					buf.str(lineA);
					buf >> gridclass.wallxsize[0];
					buf >> gridclass.wallxsize[1];
					buf >> gridclass.wallxsize[2];
					buf >> gridclass.wallxsize[3];
					buf >> gridclass.wallxsize[4];
					buf >> gridclass.wallxsize[5];
					buf.clear();
					getline(f, lineA);
					buf.str(lineA);
					buf >> gridclass.wallx_nx;
					buf >> gridclass.wallx_ny;
					buf >> gridclass.wallx_nz;
					buf.clear();
					getline(f, lineA);
					buf.str(lineA);
					buf >> gridclass.dt;
					buf.clear();
					gridclass.wallx_start[0] = gridclass.wallxsize[0];
					gridclass.wallx_start[1] = gridclass.wallxsize[2];
					gridclass.wallx_start[2] = gridclass.wallxsize[4];
					if (wavetype == 1) // if irregular wave
						wavetype = 3; // change to special case irregular with gridded wallx condition
				}
			}
		}

		if (!lineA.compare("[lsgrid]")) {
			// Lagrangian streched grid
			getline(f, lineA);
			buf.str(lineA);
			buf >> sgrid.domain[0];
			buf >> sgrid.domain[1];
			buf >> sgrid.domain[2];
			buf >> sgrid.domain[3];
			std::cout << "LS grid domain bounds: " << std::endl << lineA << std::endl;

			buf.clear();
			getline(f, lineA);
			buf.str(lineA);
			buf >> sgrid.nx;
			buf >> sgrid.ny;
			if (numparams(lineA) == 3) {
				buf >> sgrid.nl;
			}
			std::cout << "Grid resolution: " << std::endl;
			std::cout << "nx: " << sgrid.nx << std::endl;
			std::cout << "ny: " << sgrid.ny << std::endl;
			std::cout << "nl: " << sgrid.nl << std::endl;

			// allocate memory
			sgrid.allocate();

			buf.clear();
			getline(f, lineA);
			buf.str(lineA);
			buf >> sgrid.t0;
			buf >> sgrid.dt;
			std::cout << "time init: " << sgrid.t0 << std::endl;
			std::cout << "dt: " << sgrid.dt << std::endl;		

			buf.clear();
			getline(f, lineA);
			trim(lineA);
			if (!lineA.compare("stretch_params")) {
				getline(f, lineA);
				buf.str(lineA);
				buf >> sgrid.tan_a;
				buf >> sgrid.tan_b;
				std::cout << "Stretching parameters set to a=" << sgrid.tan_a << ", b=" << sgrid.tan_b << std::endl;
			}
				
			buf.clear();
			getline(f, lineA);
			trim(lineA);
			if (!lineA.compare("ignore_subdomain")) {
				getline(f, lineA);
				buf.str(lineA);
				buf >> sgrid.domain_ignore[0];
				buf >> sgrid.domain_ignore[1];
				buf >> sgrid.domain_ignore[2];
				buf >> sgrid.domain_ignore[3];
				std::cout << "The following subdomain will be ignored after intialization: " << std::endl << lineA << std::endl;
				sgrid.set_ignore();
			}
			buf.clear();

			if (wavetype == 1)
				wavetype = 4;
			else {
				std::cerr << "LS grid may only be used with irregular waves" << std::endl;
				exit(-1);
			}

		}
	}
	std::cout << "Input file read successfully." << std::endl;

	if (wavetype == 1 || wavetype == 2 || wavetype == 3 || wavetype == 4) {
		irregular.depth = depth;
		irregular.mtheta = mtheta;
		irregular.tofmax = tofmax;
		irregular.fpoint[0] = x_pos;
		irregular.fpoint[1] = y_pos;
		irregular.swl = swl;
		irregular.normalize_data();

	}
	else if (wavetype == 21) {
		stokes5.depth = depth;
		stokes5.theta = mtheta;
		stokes5.x0 = x_pos;
		stokes5.y0 = y_pos;
		stokes5.t0 = tofmax;
		stokes5.z0 = swl;
		// set the properties of the wave
		stokes5.set_stokes5_properties(wave_length, wave_height);
	}

	if (wavetype == 4) {
		sgrid.water_depth = depth;
		sgrid.initialize_surface_elevation(irregular,sgrid.t0);
		sgrid.initialize_kinematics(irregular);
	}

	return 0;

}

/*
int read_inputdata()
{
	std::string lineA;
	std::ifstream fid;
	std::string res;

	// READ INPUT FILE AND REMOVE COMMENT LINES
	fid.open("./waveinput.dat");

	// check one step up in the folder tree (this is used in the latest comflow version)
	if (fid.fail()) {
		fid.open("../waveinput.dat");
	}
	// Error check
	if (fid.fail()) {
		std::cerr << "Could not open file (is it really there?) " << std::endl;
		return -1;
		exit(1);
	}
	else {
		std::cout << "Reading data from file: waveinput.dat..." << std::endl;
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

*/



double wave_VeloX(double xpt, double ypt, double zpt, double tpt)
{

	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-depth, zpt);

	switch (wavetype) {
	// irregular waves
	case 1:
		return ramp.ramp(tpt,xpt,ypt)*irregular.u(tpt, xpt, ypt, zpt);
	// irregular gridded waves
	case 2:
		if (gridclass.initkin == 0) {
			std::cout << "Generating kinematics for interpolation:" << std::endl;
			gridclass.initialize_kinematics(irregular, 0.0);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.u(xpt, ypt, zpt);
	
		// irregular waves generated from boundary, propagate into still water
	case 3:
		// check timestep and see if updating is required
		if (!gridclass.CheckTime(tpt)) {
			gridclass.update_boundary_wallx(irregular, tpt);
		}
		// check if coordinate is within boundary wallbox
		if (!gridclass.CheckBounds(gridclass.wallxsize, xpt, ypt, zpt)) {
			gridclass.redefine_boundary_wallx(irregular, tpt, xpt, ypt, zpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.u_wall(tpt, xpt, ypt, zpt);
		
		// irregular LSgrid waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
			#pragma omp single
			sgrid.update(irregular, tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * sgrid.u(tpt, xpt, ypt, zpt);

		
	// stokes 5th
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.u(tpt, xpt, ypt, zpt);
	
	case 11:
		return ramp.ramp(tpt, xpt, ypt) * wavemaker.u_piston(tpt);
	default:
		return 0.0;
	}


}

double wave_VeloX_slope(double xpt, double ypt, double zpt, double tpt, double sx)
{

	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-depth, zpt);

	switch (wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.u(tpt, xpt, ypt, zpt);
		// irregular gridded waves
	case 2:
		if (gridclass.initkin == 0) {
			std::cout << "Generating kinematics for interpolation:" << std::endl;
			gridclass.initialize_kinematics(irregular, 0.0);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.u(xpt, ypt, zpt);

		// irregular LSgrid waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.u(tpt, xpt, ypt, zpt);
		// stokes 5th
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.u(tpt, xpt, ypt, zpt);

	case 11:
		return ramp.ramp(tpt, xpt, ypt) * wavemaker.u_piston(tpt);
	default:
		return 0.0;
	}


}


double wave_VeloY(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-depth, zpt);

	switch (wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.v(tpt, xpt, ypt, zpt);
		// irregular gridded waves
	case 2:
		if (gridclass.initkin == 0) {
			std::cout << "Generating kinematics for interpolation:" << std::endl;
			gridclass.initialize_kinematics(irregular,0.0);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.v(xpt, ypt, zpt);
	case 3:
		// check timestep and see if updating is required
		if (!gridclass.CheckTime(tpt)) {
			gridclass.update_boundary_wallx(irregular,tpt);
		}
		// check if coordinate is within boundary wallbox
		if (!gridclass.CheckBounds(gridclass.wallxsize, xpt, ypt, zpt)) {
			gridclass.redefine_boundary_wallx(irregular, tpt, xpt, ypt, zpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.v_wall(tpt, xpt, ypt, zpt);
		// irregular LSgrid waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
			#pragma omp single
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.v(tpt, xpt, ypt, zpt);
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.v(tpt, xpt, ypt, zpt);

	default:
		return 0.0;
	}
}


double wave_VeloZ(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-depth, zpt);

	switch (wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.w(tpt, xpt, ypt, zpt);
		// irregular gridded waves
	case 2:
		if (gridclass.initkin == 0) {
			std::cout << "Generating kinematics for interpolation:" << std::endl;
			gridclass.initialize_kinematics(irregular, 0.0);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.w(xpt, ypt, zpt);
	case 3:
		// check timestep and see if updating is required
		if (!gridclass.CheckTime(tpt)) {
			gridclass.update_boundary_wallx(irregular, tpt);
		}
		// check if coordinate is within boundary wallbox
		if (!gridclass.CheckBounds(gridclass.wallxsize, xpt, ypt, zpt)) {
			gridclass.redefine_boundary_wallx(irregular, tpt, xpt, ypt, zpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.w_wall(tpt, xpt, ypt, zpt);
		// irregular LSgrid waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.w(tpt, xpt, ypt, zpt);
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.w(tpt, xpt, ypt, zpt);

	default:
		return 0.0;
	}
}

// todo: not fully functional, but not really used either.... should be updated at some point
double wave_DynPres(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-depth, zpt);

	switch (wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.dp(tpt, xpt, ypt, zpt);
		// irregular gridded waves
	case 2:
		return 0.;
	case 3:
		return 0.;

	default:
		return 0.0;
	}
}

//
double wave_SurfElev(double xpt, double ypt, double tpt)
{
	switch (wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		//return waveelev(tpt, xpt, ypt);
		return ramp.ramp(tpt, xpt, ypt) * irregular.eta(tpt, xpt, ypt);
		//return irregular.eta(tpt, xpt, ypt);
		// Linear wave theory, constant profile used above free surface
	case 2:
		if (gridclass.initsurf == 0) {
			std::cout << "Initializing surface elevation storage:" << std::endl;
			gridclass.initialize_surface_elevation(irregular, 0.0);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.eta(xpt, ypt);
	case 3:
		if (!gridclass.initialized) {
			gridclass.allocate_wallx_memory();
			gridclass.init_boundary_wallx(irregular, tpt);
			gridclass.initialized = true;
		}
		// check timestep and see if updating is required
		if (!gridclass.CheckTime(tpt)) {
			gridclass.update_boundary_wallx(irregular, tpt);
			std::cout << "updating wall X kinematics tables for time interval t0=" << gridclass.t0 + gridclass.dt << "sec, t1=" << gridclass.t1 + gridclass.dt << "sec." << std::endl;
		}
		// check if coordinate is within boundary wallbox
		if (!gridclass.CheckBounds(gridclass.wallxsize, xpt, ypt, swl)) {
			gridclass.redefine_boundary_wallx(irregular, tpt, xpt, ypt, swl);
		}
		return ramp.ramp(tpt, xpt, ypt) * gridclass.eta_wall(tpt, xpt, ypt);
	case 4:
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.eta(tpt, xpt, ypt);
	case 11:
		return ramp.ramp(tpt, xpt, ypt) * wavemaker.wave_elev_piston(tpt);
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.eta(tpt, xpt, ypt);
		
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
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "CFD WAVEMAKER v.2.1.1" << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	
	// Check if license has expired
	if (check_license() == 1){
		std::cout << "License checks out...carry on..." << std::endl;
	}
	else {
		std::cout << "License for CFDwavemaker has expired. Please contact Oeystein Lande to update the program." << std::endl << std::endl;
		std::cout << "This program will auto-distruct in \n5..." << std::endl;
		wait(1);
		std::cout << "4..." << std::endl;
		wait(1);
		std::cout << "3..." << std::endl;
		wait(1);
		std::cout << "2..." << std::endl;
		wait(1);
		std::cout << "1..." << std::endl;
		wait(1);
		std::cout << "Bang!" << std::endl;

		return -1;
	}
	//for (int i = 0; i < nfreq; i++) {
	//	k[i] = pow(2. * pi * f[i], 2.) / 9.81;
	//}
	int i = read_inputdata_v2();

	return 0;
}


int wave_Cleanup()
{
	/*if (wavetype == 1) {
		delete irregular;
	}
	else if (wavetype == 4) {
		delete[] PD_time, PD_ampl, PD_velo, PD_eta;
	}*/
	return 0;
}

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

double wave_phase_velocity(int opt) {
	return irregular.phase_velocity(opt);
}

double wave_mean_length(int opt) {
	return irregular.mean_wave_length(opt);
}
double wave_mean_period(int opt) {
	return irregular.mean_wave_period(opt);
}

void wave_sgrid_update(double tpt) {
	if (!sgrid.CheckTime(tpt)) {
#pragma omp single
		sgrid.update(irregular, tpt);
	}
}

/*
int main() {
	//cout << GetCurrentWorkingDir() << endl;
	read_inputdata_v2();
	//for (int i = 0; i < 1; i++) {
	//	std::cout << double(i) << "wave elevation: " << wave_SurfElev(0.0, double(i), 0.0) << std::endl;
	//}

	std::cout << "rampvalue: " << ramp.ramp_init_y_up << std::endl;
	std::cout << "wave elevation: " << wave_SurfElev(0.0, 11.0, 0.0) << std::endl;
	std::cout << "wave elevation true: " << irregular.eta(0.0, 0.0, 0.0) << std::endl;
	std::cout << "velo x: " << wave_VeloX(0.0, 0.0, -5.0, 10.5) << std::endl;
	//std::cout << "velo y: " << wave_VeloY(0.0, 0.0, -5.0, 10.5) << std::endl;
	//std::cout << "velo z: " << wave_VeloZ(0.0, 0.0, -5.0, 10.5) << std::endl;
	//std::cout << irregular.Ampspec[0] << std::endl;

}
*/