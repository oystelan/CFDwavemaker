//
// CFDwavemaker - an Open-Source wave kinematics library
//
// This program has the soul purpose of providing wave kinematics input to any type
// of CFD program which can be linked up as a dynamic link library or statically.
// The program is created and updated by the Oeystein Lande. The link library
// compiles on both linux and windows. The following wave theory/types are currently
// supported:
// linear wave theory (2D and 3D), second order wave theory (2D and 3D) (sharma &
// dean), wave paddle theory (2D only at the moment)
//
//
// Current version: v214
// Date: 2020-12-05
// (c) Oystein Lande
// --------------------------------------------------------------------------------
#include <stdio.h>
//#include <cstdlib>
//#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <limits>
#include <ctime>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
//#include <filesystem>


//#include <direct.h> // windows only function
#//include <cctype>
//#include <locale>
#include "CFDwavemaker.h" 
#include "Stokes5.h"
#include "Irregular.h"
#include "Utils.h"
#include "Wavemaker.h"
#include "lsgrid.h"
#include "lsgrid_spline.h"
#include "probes.h"

// for now, spectralwavedata is only supported in linux build

#if defined(SWD_enable)
#include "SpectralWaveData.h"
#endif

#if defined(VTK_enable)
#include "VTKreader.h"
#endif

#define largeval 1.E12
//#include <fftw3.h>


// Variables
//int nfreq, ndir, wavetype, extmet, pertmet, meth, bandwidth, n_timesteps, rampswitch, normalizeA, spreadfunc;




// Storage class of input data.
class CFDwavemakerInputdata {
public:
	double depth;
	double x_pos, y_pos, tofmax, current_speed, wave_length, wave_height;
	double mtheta = 0.;
	double swl = 0.;
	bool bw_auto_calc = false;
	int wavetype;
	bool property_read = false;
	double rho = 1025.;
	double gravity = 9.81;
	// SWD parameters
	double nsumx = -1, nsumy = -1, impl = 0, ipol = 0, norder = 0;
	bool dc_bias = false;
	std::string swdFileName;

	// lsgrid related data
	int lsgrid_interp_scheme = 0; // 0=linear, 1=spline
	double lsgrid_domain[4] = {};
	int lsgrid_nx = 5;
	int lsgrid_ny = 0;
	int lsgrid_nl = 4;
	double lsgrid_t0 = 0.;
	double lsgrid_dt = 0.1;
	double lsgrid_tan_a;
	double lsgrid_tan_b;
	double lsgrid_domain_ignore[4] = {largeval, -largeval, largeval, -largeval};
	bool lsgrid_ignore_domain = false;
	bool lsgrid_ignore_at_init = 0;
	bool lsgrid_init_only = 0;
	std::string lsgrid_scheme_str;

	CFDwavemakerInputdata() {
	}

	~CFDwavemakerInputdata() {
	}
};

bool CFDwmInit = false;

//double ampl, depth, s, mtheta, tofmax, fpoint[2], trampdata[3], xrampdata[3], yrampdata[3];

CFDwavemakerInputdata inputdata;

// Stokes 5 class
Stokes5 stokes5;

// Irregular class
Irregular irregular;

// Wavemaker theory class
Wavemaker wavemaker;

// Grid class
lsGrid sgrid;
lsGridSpline sgrids;

// Ramp class
Ramp ramp;

// probes class
Probes probes;

// SWD class;
#if defined(SWD_enable)
SpectralWaveData *swd;
#endif

#if defined(VTK_enable)
VTKreader vtkreader;
#endif

//string GetCurrentWorkingDir(void) {
//	char buff[FILENAME_MAX];
//	GetCurrentDir(buff, FILENAME_MAX);
//	std::string current_working_dir(buff);
//	return current_working_dir;
//}

// Some useful utilitize functions

void wait(int seconds)
{
	clock_t endwait;
	endwait = clock() + seconds * CLOCKS_PER_SEC;
	while (clock() < endwait) {}
}

/* A sorting function for vectors whihc returns the indices after sorting */
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
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

/* main input file reader function*/
int process_inputdata(std::string res, Irregular& irreg, Stokes5& stokes, Wavemaker& wmaker, lsGrid& lsgrid, lsGridSpline& lsgrids, Ramp& rramp) {
	std::string lineA, dummystr, lineP;
	std::ifstream fid;
	std::istringstream buf;
	std::istringstream f(res);
	bool skip_getline = false;

	std::cout << std::boolalpha; // display booleans as true or false when printed (instead of 0 or 1)

	//get and write data lines
	while (!f.eof()) {
		if (skip_getline) {
			skip_getline = false;
		}
		else {
			getline(f, lineA);
			trim(lineA);
		}
		
		//std::cout << lineA << std::endl;
		// Convension for internal wavetype numbering:
		// irregular wave theory variants: 1-10
		// wavemaker theory variants: 11-20
		// Regular wave theories: 21-30
		// HOSM and other: 31-40
		// VTK: 41-50
		if (!lineA.compare("[wave type]")) {
			std::cout << "----------" << std::endl;
			std::cout << "Wave type:" << std::endl;
			std::cout << "----------" << std::endl;
			getline(f, lineA);
			trim(lineA);
			//std::cout << lineA << std::endl;
			// check if valid wave type is given
			if (!lineA.compare("irregular") || !lineA.compare("1")) {
				inputdata.wavetype = 1;
				std::cout << "Irregular perturbation wave theory specified" << std::endl;
			}
			else if (!lineA.compare("wavemaker")) {
				inputdata.wavetype = 11;
				std::cout << "Wave maker theory specified" << std::endl;
			}
			else if (!lineA.compare("regular")) {
				inputdata.wavetype = 21;
				std::cout << "Regular 5th order Stokes wave specified" << std::endl;
			}
			else if (!lineA.compare("swd")) {
				inputdata.wavetype = 31;
				std::cout << "Spectral wave data (swd) specified" << std::endl;
			}
			else if (!lineA.compare("vtk")) {
				inputdata.wavetype = 41;
				std::cout << "Special branch VTK file interpolation" << std::endl;
			}
			
			else {
				std::cerr << "INPUTFILE ERROR: Unknown wave type specified. Valid alternatives are:" << std::endl;
				std::cerr << "irregular" << std::endl;
				std::cerr << "regular" << std::endl;
				std::cerr << "wavemaker" << std::endl;
				std::cerr << "swd" << std::endl;
				std::cerr << "vtk" << std::endl;
				
				exit(-1);
			}
		}

		if (!lineA.compare("[general input data]")) { //mandatory
			std::cout << "-------------------" << std::endl;
			std::cout << "General input data:" << std::endl;
			std::cout << "-------------------" << std::endl;
			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 5, "depth")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.depth;
					buf.clear();
					std::cout << "Water depth: " << inputdata.depth << "m" << std::endl;
				}
				if (!lineA.compare(0, 6, "mtheta")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.mtheta;
					buf.clear();
					std::cout << "Mean wave direction: " << inputdata.mtheta << "degrees" << std::endl;
				}
				if (!lineA.compare(0, 9, "normalize")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> irreg.normalize;
					buf.clear();
					std::cout << "Normalize wave amplitudes by spectral zeroth moment: " << irreg.normalize << std::endl;
				}
				if (!lineA.compare(0, 7, "amplify")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> irreg.ampl;
					buf.clear();
					std::cout << "Amplify (gain): " << irreg.ampl << "m" << std::endl;
				}
				if (!lineA.compare(0, 3, "swl")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.swl;
					buf.clear();
					std::cout << "Still water line set to: " << inputdata.swl << "m" << std::endl;
				}
				if (!lineA.compare(0, 7, "gravity")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.gravity;
					buf.clear();
					std::cout << "Gravity set to: " << inputdata.gravity << "m/s^2" << std::endl;
				}
				if (!lineA.compare(0, 3, "rho")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.rho;
					buf.clear();
					std::cout << "Water density rho set to: " << inputdata.gravity << "kg/m^3" << std::endl;
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
			
		}
		
		if (!lineA.compare("[second order]")) { //optional
			
			std::cout << "Second order irregular waves Specified." << std::endl;

			if (inputdata.wavetype > 10) {
				std::cerr << "INPUTFILE ERROR: This tag is only available for irregular waves." << std::endl;
				exit(-1);
			}		
			irreg.order = 2;

			while (!f.eof()) {
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 6, "extval")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> irreg.extrapolation_met;
					buf.clear();
					std::cout << "Second order velocity extrapolation method for z > swl: " << irreg.extrapolation_met << std::endl;
				}
				if (!lineA.compare(0, 6, "cutoff")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> irreg.dw_cutoff;
					buf.clear();
					std::cout << "Second order high-frequency cut-off set to: " << irreg.dw_cutoff << std::endl;
				}
				if (!lineA.compare(0, 9, "bandwidth")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;
					
					if (!dummystr.compare(0, 3, "off")) {
						// Do nothing. default value is already a very high number
						std::cout << "Bandwidth: off" << std::endl;
					}
					else if (!dummystr.compare(0, 4, "auto")) {
						// Compute a decent bandwidth value. todo: make a function which does this
						std::cout << "Bandwidth: auto" << std::endl;
						inputdata.bw_auto_calc = true;
					}
					else { // assumes that a value is given
						irreg.dw_bandwidth = atof(dummystr.c_str());
						std::cout << "Bandwidth: " << irreg.dw_bandwidth << " rad/s" << std::endl;
					}

					buf.clear();
				}			
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}

		if (!lineA.compare("[wave reference point]")) { //optional
			std::cout << "---------------------" << std::endl;
			std::cout << "Wave reference point:" << std::endl;
			std::cout << "---------------------" << std::endl;
			while (!f.eof()) {
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 4, "time")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.tofmax;
					buf.clear();
					std::cout << "t0: " << inputdata.tofmax << " sec" <<std::endl;
				}
				if (!lineA.compare(0, 1, "x")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.x_pos;
					buf.clear();
					std::cout << "x0: " << inputdata.x_pos << " m" << std::endl;
				}
				if (!lineA.compare(0, 1, "y")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.y_pos;
					buf.clear();
					std::cout << "y0: " << inputdata.y_pos << " m" << std::endl;
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}

		if (!lineA.compare("[ramps]")) { //optional
			std::cout << "------" << std::endl;
			std::cout << "Ramps:" << std::endl;
			std::cout << "------" << std::endl;
			rramp.ramp_init = true;
			// read time ramp data
			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 11, "time_rampup")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> rramp.ramp_init_time_up;
					buf >> rramp.time_rampup_start;
					buf >> rramp.time_rampup_end;
					buf.clear();
					std::cout << "time-rampup: " << rramp.ramp_init_time_up << ", starttime: " << rramp.time_rampup_start  << " sec, endtime: " << rramp.time_rampup_end << " sec." << std::endl;
				}
				if (!lineA.compare(0, 13, "time_rampdown")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> rramp.ramp_init_time_down;
					buf >> rramp.time_rampdown_start;
					buf >> rramp.time_rampdown_end;
					buf.clear();
					std::cout << "time-rampdown: " << rramp.ramp_init_time_down << ", starttime: " << rramp.time_rampdown_start << " sec, endtime: " << rramp.time_rampdown_end << " sec." << std::endl;
				}
				if (!lineA.compare(0, 8, "x_rampup")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> rramp.ramp_init_x_up;
					buf >> rramp.x_rampup_start;
					buf >> rramp.x_rampup_end;
					buf.clear();
					std::cout << "x-rampup: " << rramp.ramp_init_x_up << ", startpos: " << rramp.x_rampup_start << " m, endpos: " << rramp.x_rampup_end << " m." << std::endl;
				}
				if (!lineA.compare(0, 10, "x_rampdown")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> rramp.ramp_init_x_down;
					buf >> rramp.x_rampdown_start;
					buf >> rramp.x_rampdown_end;
					buf.clear();
					std::cout << "x-rampdown: " << rramp.ramp_init_x_down << ", startpos: " << rramp.x_rampdown_start << " m, endpos: " << rramp.x_rampdown_end << " m." << std::endl;
				}
				if (!lineA.compare(0, 8, "y_rampup")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> rramp.ramp_init_y_up;
					buf >> rramp.y_rampup_start;
					buf >> rramp.y_rampup_end;
					buf.clear();
					std::cout << "y-rampup: " << rramp.ramp_init_y_up << ", startpos: " << rramp.y_rampup_start << " m, endpos: " << rramp.y_rampup_end << " m." << std::endl;
				}
				if (!lineA.compare(0, 10, "y_rampdown")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> rramp.ramp_init_y_down;
					buf >> rramp.y_rampdown_start;
					buf >> rramp.y_rampdown_end;
					buf.clear();
					std::cout << "y-rampdown: " << rramp.ramp_init_y_down << ", startpos: " << rramp.y_rampdown_start << " m, endpos: " << rramp.y_rampdown_end << " m." << std::endl;
				}
				
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}
		// Wave properties: irregular wave components (manual)
		if (!lineA.compare("[irregular wave components]")) {
			std::cout << "--------------------------" << std::endl;
			std::cout << "Irregular Wave components:" << std::endl;
			std::cout << "--------------------------" << std::endl;

			if (inputdata.wavetype != 1) {
				std::cerr << "INPUTFILE ERROR: irregular wave components does not match the specified wave type. Check inputfile" << std::endl;
				exit(-1);
			}
			// if (irreg.initialized) {
			// 	std::cerr << "INPUTFILE ERROR: Irregular wave already initialized. please check input file" << std::endl;
			// 	exit(-1);
			// }
			getline(f, lineA);
			trim(lineA);
			if (!lineA.compare(0, 5, "nfreq")) {
				buf.str(lineA);
				buf >> dummystr;
				buf >> irreg.nfreq;
				buf.clear();
			}
			else {
				std::cerr << "INPUTFILE ERROR: parameter nfreq missing or not specified correctly" << std::endl;
				exit(-1);
			}
			getline(f, lineA);
			trim(lineA);
			if (!lineA.compare(0, 4, "ndir")) {
				buf.str(lineA);
				buf >> dummystr;
				buf >> irreg.ndir;
				buf.clear();
			}
			else {
				std::cerr << "INPUTFILE ERROR: parameter ndir missing or not specified correctly" << std::endl;
				exit(-1);
			}
			// if ndir = 0. a single direction is read for each frequency.
			if (irreg.ndir == 0) {
				std::cout << "Irregular seas, one directional component for each frequency specified" << std::endl;
				std::cout << "Number of frequency components: " << irreg.nfreq << std::endl;
				irreg.ndir = 1;

				// Create some temporary vectors for storage of spectral data
				std::vector<double> omega;
				std::vector<double> A;
				std::vector<double> k;
				std::vector<double> theta;
				std::vector<double> phase;

				double temp;
				std::cout << "# OMEGA[rad / s]    A[m]           K             Phase[rad]        Theta [rad]" << std::endl;
				for (int i = 0; i < irreg.nfreq; i++) {
					getline(f, lineA);
					// if new tag is reach. break while loop.
					if (!lineA.compare(0, 1, "[")) {
						std::cerr << "INPUTFILE ERROR: Not enough frequency components specified. please check to make sure that the number of lines matches the specified number of frequency components" << std::endl;
						exit(-1);
					}
					std::cout << lineA << std::endl;
					buf.str(lineA);
					buf >> temp;
					omega.push_back(temp);
					buf >> temp;
					A.push_back(temp);
					buf >> temp;
					k.push_back(temp);
					buf >> temp;
					phase.push_back(temp);
					buf >> temp;
					theta.push_back(temp);
					buf.clear();
				}

				// Sort vectors as a function of omega (ascending)
				for (auto i : sort_indices(omega)) {
					//std::cout << omega[i] << std::endl;
					irreg.omega.push_back(omega[i]);
					irreg.A.push_back(A[i]);
					irreg.k.push_back(k[i]);
					irreg.phase.push_back(phase[i]);
					irreg.theta.push_back(theta[i]);
				}
				irreg.initialized = true;

			}
			// A spreading function is used
			else {
				std::cout << "Irregular seas, directional spreading defined separately" << std::endl;
				std::cout << "Number of frequency components: " << irreg.nfreq << std::endl;
				std::cout << "Number of directional components: " << irreg.ndir << std::endl;
				// Read frequency data (omega, Sw and K)
				double* omega_temp = new double[irreg.nfreq];
				double* Ampspec_temp = new double[irreg.nfreq];
				double* k_temp = new double[irreg.nfreq];
				double* phas_temp = new double[irreg.nfreq];
				std::cout << "# OMEGA[rad / s]    A[m]           K             Phase[rad]" << std::endl;
				for (int i = 0; i < irreg.nfreq; i++) {
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
				double* theta_temp = new double[irreg.ndir];
				double* D_ampl_temp = new double[irreg.ndir];
				std::cout << "# Theta [rad] D(theta)" << std::endl;
				for (int i = 0; i < irreg.ndir; i++) {
					getline(f, lineA);
					if (!lineA.compare(0, 1, "[")) {
						std::cerr << "INPUTFILE ERROR: too few directional components.." << std::endl;
						exit(-1);
					}
					buf.str(lineA);
					buf >> theta_temp[i];
					buf >> D_ampl_temp[i];
					buf.clear();
				}

				// Create some temporary vectors for storage of spectral data
				std::vector<double> omega;
				std::vector<double> A;
				std::vector<double> k;
				std::vector<double> theta;
				std::vector<double> phase;

				for (int i = 0; i < irreg.nfreq; i++) {
					for (int j = 0; j < irreg.ndir; j++) {
						omega.push_back(omega_temp[i]);
						k.push_back(k_temp[i]);
						A.push_back(Ampspec_temp[i] * D_ampl_temp[j]);
						phase.push_back(phas_temp[i]);
						theta.push_back(theta_temp[j]);

					}
				}
				delete[] Ampspec_temp, omega_temp, phas_temp, k_temp, theta_temp, D_ampl_temp;

				// Sort vectors as a function of omega (ascending)
				for (auto i : sort_indices(omega)) {
					//std::cout << omega[i] << std::endl;
					irreg.omega.push_back(omega[i]);
					irreg.A.push_back(A[i]);
					irreg.k.push_back(k[i]);
					irreg.phase.push_back(phase[i]);
					irreg.theta.push_back(theta[i]);
				}
				irreg.initialized = true;
			}
			inputdata.property_read = true;
		}

		// Wave properties: piston wave maker
		if (!lineA.compare("[pistonwavemaker wave properties]")) {
			if (inputdata.wavetype != 11) {
				std::cerr << "piston wave maker properties does not match the specified wave type. Check inputfile" << std::endl;
				exit(-1);
			}
			if (wmaker.initialized) {
				std::cerr << "INPUTFILE ERROR: Pistonwavemaker already initialized. please check input file" << std::endl;
				exit(-1);
			}
			// Wavemaker theory (piston)
			// read alpha values
			getline(f, lineA);
			buf.str(lineA);
			buf >> wmaker.alpha_z;
			buf >> wmaker.alpha_u;
			buf.clear();

			getline(f, lineA);
			buf.str(lineA);
			buf >> wmaker.n_timesteps;
			//n_timesteps = stoi(lineA);
			std::cout << "Number of timesteps: " << inputdata.wavetype << std::endl;

			// declare some vectors to store piston data
			wmaker.PD_time = new double[wmaker.n_timesteps];
			wmaker.PD_ampl = new double[wmaker.n_timesteps];
			wmaker.PD_velo = new double[wmaker.n_timesteps];
			wmaker.PD_eta = new double[wmaker.n_timesteps];

			for (int i = 0; i < wmaker.n_timesteps; i++) {
				getline(f, lineA);
				buf.str(lineA);
				buf >> wmaker.PD_time[i];
				buf >> wmaker.PD_ampl[i];
				buf >> wmaker.PD_velo[i];
				buf >> wmaker.PD_eta[i];
				buf.clear();
			}

			inputdata.property_read = true;
			wmaker.initialized = true;
		}
		
		// Wave properties: Stokes wave properties
		if (!lineA.compare("[stokes wave properties]")) {
			if (inputdata.wavetype != 21) {
				std::cerr << "INPUTFILE ERROR: Stokes wave properties does not match the specified wave type. Check inputfile" << std::endl;
				exit(-1);
			}
			if (stokes.initialized) {
				std::cerr << "INPUTFILE ERROR: Stokes 5th wave already initialized. please check input file" << std::endl;
				exit(-1);
			}
			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 11, "wave_length")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> stokes.wave_length;
					buf.clear();
					std::cout << "Wave length:  " << stokes.wave_length << std::endl;
				}
				if (!lineA.compare(0, 11, "wave_height")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> stokes.wave_height;
					buf.clear();
					std::cout << "Wave height:  " << stokes.wave_height << std::endl;
				}
				if (!lineA.compare(0, 13, "current_speed")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> stokes.current;
					buf.clear();
					std::cout << "current speed:  " << stokes.current << std::endl;
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
			inputdata.property_read = true;
			stokes.initialized = true;
		}

		if (!lineA.compare("[swd wave properties]")) { //optional
			std::cout << "-------------------------------------" << std::endl;
			std::cout << "Spectral wave data:" << std::endl;
			std::cout << "-------------------------------------" << std::endl;
			if (inputdata.wavetype != 31) {
				std::cerr << "INPUTFILE ERROR: please set [wave type] = swd when using keyword [swd wave properties]." << std::endl;
				exit(-1);
			}
			
			while (!f.eof()) {
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 7, "swdfile")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.swdFileName;
					buf.clear();
					std::cout << "Swd file specified:  " << inputdata.swdFileName << std::endl;
				}
				if (!lineA.compare(0, 3, "rho")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.rho;
					buf.clear();
					std::cout << "rho:  " << inputdata.rho << " kg/m^3." << std::endl;
				}
				if (!lineA.compare(0, 5, "nsumx")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.nsumx;
					buf.clear();
					std::cout << "nsumx:  " << inputdata.nsumx << std::endl;
				}
				if (!lineA.compare(0, 5, "nsumy")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.nsumy;
					buf.clear();
					std::cout << "nsumy:  " << inputdata.nsumy << std::endl;
				}
				if (!lineA.compare(0, 4, "impl")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.impl;
					buf.clear();
					std::cout << "impl:  " << inputdata.impl << std::endl;
				}
				if (!lineA.compare(0, 4, "ipol")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.ipol;
					buf.clear();
					std::cout << "ipol:  " << inputdata.ipol << std::endl;
				}
				if (!lineA.compare(0, 6, "norder")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.norder;
					buf.clear();
					std::cout << "norder:  " << inputdata.norder << std::endl;
				}
				if (!lineA.compare(0, 9, "dc_bias")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;

					if (!dummystr.compare(0, 5, "false")) {
						// Do nothing. default value is already a very high number
						inputdata.dc_bias = false;
					}
					else if (!dummystr.compare(0, 4, "true")) {
						// Compute a decent bandwidth value. todo: make a function which does this						
						inputdata.dc_bias = true;
					}
					else { // assumes that a value is given
						inputdata.dc_bias = atof(dummystr.c_str());
					}
					std::cout << "dc_bias:   " << inputdata.dc_bias << std::endl;

					buf.clear();
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}		
		}


		if (!lineA.compare("[lsgrid]")) {
			std::cout << "-----------------------------------" << std::endl;
			std::cout << "Lagrangian Stretched grid (lsgrid):" << std::endl;
			std::cout << "-----------------------------------" << std::endl;			
			
			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				std::cout << lineA << std::endl;
				trim(lineA);
				if (!lineA.compare(0, 6, "scheme")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;
					inputdata.lsgrid_scheme_str = dummystr;
					if (!dummystr.compare(0, 6, "linear")) {
						inputdata.lsgrid_interp_scheme = 0;
						// Do nothing. default value is already a very high number
						std::cout << "Linear interpolation scheme (default) chosen." << std::endl;
					}
					else if (!dummystr.compare(0, 6, "spline")) {
						inputdata.lsgrid_interp_scheme = 1;
						// Compute a decent bandwidth value. todo: make a function which does this
						std::cout << "Spline interpolation scheme chosen." << std::endl;
					}
					else { // assumes that a value is given
						std::cout << "Unknown interpolation scheme specified (use either linear or spline)." << std::endl;
						exit(-1);
					}

					buf.clear();
				}		
				if (!lineA.compare(0, 6, "bounds")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_domain[0];
					buf >> inputdata.lsgrid_domain[1];
					buf >> inputdata.lsgrid_domain[2];
					buf >> inputdata.lsgrid_domain[3];
					buf.clear();
				}
				if (!lineA.compare(0, 2, "nx")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_nx;
					buf.clear();
				}
				if (!lineA.compare(0, 2, "ny")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_ny;
					buf.clear();
				}
				if (!lineA.compare(0, 2, "nl")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_nl;
					buf.clear();
				}
				if (!lineA.compare(0, 2, "t0")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_t0;
					buf.clear();
				}
				if (!lineA.compare(0, 2, "dt")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_dt;
					buf.clear();
				}
				if (!lineA.compare(0, 14, "stretch_params")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_tan_a;
					buf >> inputdata.lsgrid_tan_b;
					buf.clear();
				}
				if (!lineA.compare(0, 16, "ignore_subdomain")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_domain_ignore[0];
					buf >> inputdata.lsgrid_domain_ignore[1];
					buf >> inputdata.lsgrid_domain_ignore[2];
					buf >> inputdata.lsgrid_domain_ignore[3];
					buf.clear();
					inputdata.lsgrid_ignore_domain = true;
				}
				if (!lineA.compare(0, 14, "ignore_at_init")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_ignore_at_init;
					buf.clear();
				}
				if (!lineA.compare(0, 9, "init_only")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.lsgrid_init_only;
					buf.clear();
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
				
			}

			// allocate memory
			

			std::cout << "Interpolation Scheme: " << inputdata.lsgrid_scheme_str << std::endl;
			std::cout << "LS grid domain bounds: " << std::endl;
			std::cout << "\tXMIN: " << inputdata.lsgrid_domain[0] << std::endl;
			std::cout << "\tXMAX: " << inputdata.lsgrid_domain[1] << std::endl;
			std::cout << "\tYMIN: " << inputdata.lsgrid_domain[2] << std::endl;
			std::cout << "\tYMAX: " << inputdata.lsgrid_domain[3] << std::endl;
			std::cout << "Grid resolution: " << std::endl;
			std::cout << "nx: " << inputdata.lsgrid_nx << std::endl;
			std::cout << "ny: " << inputdata.lsgrid_ny << std::endl;
			std::cout << "nl: " << inputdata.lsgrid_nl << std::endl;
			std::cout << "time init: " << inputdata.lsgrid_t0 << std::endl;
			std::cout << "dt: " << inputdata.lsgrid_dt << std::endl;
			std::cout << "Stretching parameters set to a=" << inputdata.lsgrid_tan_a << ", b=" << inputdata.lsgrid_tan_b << std::endl;
			if (inputdata.lsgrid_ignore_domain) {
				std::cout << "The following subdomain will be ignored after intialization: " << std::endl;
				std::cout << "\tXMIN: " << inputdata.lsgrid_domain_ignore[0] << std::endl;
				std::cout << "\tXMAX: " << inputdata.lsgrid_domain_ignore[1] << std::endl;
				std::cout << "\tYMIN: " << inputdata.lsgrid_domain_ignore[2] << std::endl;
				std::cout << "\tYMAX: " << inputdata.lsgrid_domain_ignore[3] << std::endl;
			}
			std::cout << "Ignore subdomain at t_init: " << inputdata.lsgrid_ignore_at_init << std::endl;

			if (inputdata.lsgrid_interp_scheme == 0) {
				if (inputdata.wavetype == 1){
					inputdata.wavetype = 4;
					std::cout << "LS grid + linear interpolation + irregular" << std::endl;
				}
				else if (inputdata.wavetype == 31){
					inputdata.wavetype = 34;
				std::cout << "LS grid + swd = True." << std::endl;
				}
				std::cout << "LS grid + irregular = True." << std::endl;
				// copy necessary parameters to lsgrid class
				lsgrid.domain[0] = inputdata.lsgrid_domain[0];
				lsgrid.domain[1] = inputdata.lsgrid_domain[1];
				lsgrid.domain[2] = inputdata.lsgrid_domain[2];
				lsgrid.domain[3] = inputdata.lsgrid_domain[3];
				lsgrid.nx = inputdata.lsgrid_nx;
				lsgrid.ny = inputdata.lsgrid_ny;
				lsgrid.nl = inputdata.lsgrid_nl;
				lsgrid.t0 = inputdata.lsgrid_t0;
				lsgrid.dt = inputdata.lsgrid_dt;
				lsgrid.tan_a = inputdata.lsgrid_tan_a;
				lsgrid.tan_b = inputdata.lsgrid_tan_b;
				lsgrid.domain_ignore[0] = inputdata.lsgrid_domain_ignore[0];
				lsgrid.domain_ignore[1] = inputdata.lsgrid_domain_ignore[1];
				lsgrid.domain_ignore[2] = inputdata.lsgrid_domain_ignore[2];
				lsgrid.domain_ignore[3] = inputdata.lsgrid_domain_ignore[3];
				lsgrid.ignore_domain = inputdata.lsgrid_ignore_domain;
				lsgrid.init_only = inputdata.lsgrid_init_only;
				lsgrid.ignore_at_init = inputdata.lsgrid_ignore_at_init;
				lsgrid.allocate();
			}
			else if (inputdata.lsgrid_interp_scheme == 1){
				inputdata.wavetype = 5;
				std::cout << "LS grid + Spline + irregular" << std::endl;
				// copy necessary parameters to lsgrid class
				lsgrids.domain[0] = inputdata.lsgrid_domain[0];
				lsgrids.domain[1] = inputdata.lsgrid_domain[1];
				lsgrids.domain[2] = inputdata.lsgrid_domain[2];
				lsgrids.domain[3] = inputdata.lsgrid_domain[3];
				lsgrids.nx = inputdata.lsgrid_nx;
				lsgrids.ny = inputdata.lsgrid_ny;
				lsgrids.nl = inputdata.lsgrid_nl;
				lsgrids.t0 = inputdata.lsgrid_t0;
				lsgrids.dt = inputdata.lsgrid_dt;
				lsgrids.tan_a = inputdata.lsgrid_tan_a;
				lsgrids.tan_b = inputdata.lsgrid_tan_b;
				lsgrids.domain_ignore[0] = inputdata.lsgrid_domain_ignore[0];
				lsgrids.domain_ignore[1] = inputdata.lsgrid_domain_ignore[1];
				lsgrids.domain_ignore[2] = inputdata.lsgrid_domain_ignore[2];
				lsgrids.domain_ignore[3] = inputdata.lsgrid_domain_ignore[3];
				lsgrids.ignore_domain = inputdata.lsgrid_ignore_domain;
				lsgrids.init_only = inputdata.lsgrid_init_only;
				lsgrids.ignore_at_init = inputdata.lsgrid_ignore_at_init;
				lsgrids.allocate();
			}
			else {
				std::cerr << "INPUTFILE ERROR: LS grid may currently only be used with irregular second order wave theory and spectral wave method" << std::endl;
				exit(-1);
			}
		}
		
		if (!lineA.compare("[vtk output]")) { //optional
			std::cout << "--------------------" << std::endl;
			std::cout << "VTK output settings:" << std::endl;
			std::cout << "--------------------" << std::endl;
			if (inputdata.wavetype != 4 && inputdata.wavetype != 34 && inputdata.wavetype != 5)  {
				std::cout << "InputError: VTK output only works in combination with LSgrid at the moment." << std::endl;
				exit(-1);
			}
			sgrid.dump_vtk = true;
			sgrids.dump_vtk = true;

			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 12, "storage_path")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> lsgrid.vtk_directory_path;
					lsgrids.vtk_directory_path = lsgrid.vtk_directory_path;
					buf.clear();
					std::cout << "Directory for storage of vtk files: " << lsgrid.vtk_directory_path << std::endl;
				}
				if (!lineA.compare(0, 8, "filename")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> lsgrid.vtk_prefix;
					lsgrids.vtk_prefix = lsgrid.vtk_prefix;
					buf.clear();
					std::cout << "filename prefix: " << lsgrid.vtk_prefix << std::endl;
				}
				if (!lineA.compare(0, 13, "vtk_timelabel")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> lsgrid.vtk_timelabel;
					lsgrids.vtk_timelabel = lsgrid.vtk_timelabel;
					buf.clear();
					std::cout << "time label to use in vtk specified as: " << lsgrid.vtk_timelabel << std::endl;
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}

		if (!lineA.compare("[timeseries output]")) { //optional
			std::cout << "----------------------------" << std::endl;
			std::cout << "Time-series output settings:" << std::endl;
			std::cout << "----------------------------" << std::endl;

			/*[timeseries output]
			storage_path ./ts/ 
				filename tsfile
				npos 3
				# x y z
				0.0 0.0 0.0
				3.3 5.4 - 10.
				0.0 5.4 - 10.*/

			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 12, "storage_path")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> probes.ts_path;
					buf.clear();
					std::cout << "Directory for storage of timeseries files: " << probes.ts_path << std::endl;
				}
				if (!lineA.compare(0, 8, "filename")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> probes.ts_filename;
					buf.clear();
					std::cout << "filename prefix: " << probes.ts_filename << std::endl;
				}
				if (!lineA.compare(0, 2, "t0")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> probes.t0;
					buf.clear();
				}
				if (!lineA.compare(0, 2, "dt")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> probes.dt;
					buf.clear();
				}
				if (!lineA.compare(0, 4, "npos")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> probes.ts_nprobes;
					buf.clear();
					std::cout << "Number of probes: " << probes.ts_nprobes << std::endl;
					double tempx, tempy, tempz;
					for (int i = 0; i < probes.ts_nprobes; i++) {
						getline(f, lineA);
						// if new tag is reach. break while loop.
						if (!lineA.compare(0, 1, "[")) {
							std::cerr << "INPUTFILE ERROR in [timeseries output]: Not enough timeseries coordnates specified" << std::endl;
							exit(-1);
						}
						std::cout << lineA << std::endl;
						buf.str(lineA);
						buf >> tempx;
						buf >> tempy;
						buf >> tempz;
						probes.coords.push_back({tempx, tempy, tempz});						
						buf.clear();

					}
					probes.init_probes();
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}

		if (!lineA.compare("[vtk input]")) {

#if defined(VTK_enable)
			std::cout << "----------------------------" << std::endl;
			std::cout << "VTK interpolation from file:" << std::endl;
			std::cout << "----------------------------" << std::endl;

			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 12, "storage_path")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> vtkreader.vtkfilepath;
					buf.clear();
					std::cout << "Directory where vtk files are stored: " << vtkreader.vtkfilepath << std::endl;
				}
				if (!lineA.compare(0, 8, "filename")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> vtkreader.vtk_prefix;
					buf.clear();
					std::cout << "filename prefix: " << vtkreader.vtk_prefix << std::endl;
				}
				if (!lineA.compare(0, 19, "name_velocity_field")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> vtkreader.Uname;
					buf.clear();
					std::cout << "Name of velocity scalar field: " << vtkreader.Uname << std::endl;
				}
				if (!lineA.compare(0, 7, "t_start")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> vtkreader.t_start;
					buf.clear();
					std::cout << "User specified start time (t0): " << vtkreader.Uname << std::endl;
				}
				if (!lineA.compare(0, 37, "update_height_function_every_timestep")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;

					if (!dummystr.compare(0, 5, "false")) {
						// Do nothing. default value is already a very high number
						vtkreader.recompute_betah_every_timstep = false;
					}
					else if (!dummystr.compare(0, 4, "true")) {
						// Compute a decent bandwidth value. todo: make a function which does this						
						vtkreader.recompute_betah_every_timstep = true;
					}
					else { // assumes that a value is given
						vtkreader.recompute_betah_every_timstep = atof(dummystr.c_str());
					}
					std::cout << "update_height_function_every_timestep:   " << vtkreader.recompute_betah_every_timstep << std::endl;
					buf.clear();
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
#else
			std::cerr << "VTK interpolation not supported in current compiled version of CFDwavemaker. Recompile width flag -DVTK_enable=1." << std::endl;
			exit(-1);
#endif
			
		}
	}

	// Do some checks of the specified input parameters. Some combinations are prohibited

	if (inputdata.wavetype == 4 && irregular.order != 2) {
		std::cerr << "LSgrid interpolation uses strictly second order theory. Doesnt make sence to use grid interpolation for linear theory since this is very fast anyway. Either add the tag [second order] to the input file, or remore [lsgrid]." << std::endl;
		exit(-1);

	}


	std::cout << "\n-----------------------------------------------" << std::endl;
	std::cout << "Input file read successfully." << std::endl;
	std::cout << "***********************************************\n\n" << std::endl;

	std::cout << "Wavetype: " << inputdata.wavetype << std::endl;


	if (inputdata.wavetype == 1 || inputdata.wavetype == 2 || inputdata.wavetype == 3 || inputdata.wavetype == 4) {
		irregular.depth = inputdata.depth;
		irregular.mtheta = inputdata.mtheta;
		irregular.tofmax = inputdata.tofmax;
		irregular.fpoint[0] = inputdata.x_pos;
		irregular.fpoint[1] = inputdata.y_pos;
		irregular.swl = inputdata.swl;
		irregular.g = inputdata.gravity;
		irregular.rho = inputdata.rho;

		if (inputdata.bw_auto_calc) {
			irreg.dw_bandwidth = irregular.bandwidth_estimator();
			std::cout << "BW_autocalc: Bandwidth parameter set to " << irreg.dw_bandwidth << " rad/s." << std::endl;		
		}
		irregular.normalize_data();
		irregular.calculate_bwindices();
		irregular.dumpSpectralComponents();

		if (!inputdata.property_read) {
			std::cerr << "INPUTFILE ERROR: Irregular wave selected, but no wave components/wave properties found in input file." << std::endl;
			exit(-1);
		}
	}
	else if (inputdata.wavetype == 11) {
		if (!inputdata.property_read) {
			std::cerr << "INPUTFILE ERROR: Wave maker selected, but no wavemaker property data found in input file." << std::endl;
			exit(-1);
		}
	}
	else if (inputdata.wavetype == 21) {
		stokes.depth = inputdata.depth;
		stokes.theta = inputdata.mtheta;
		stokes.x0 = inputdata.x_pos;
		stokes.y0 = inputdata.y_pos;
		stokes.t0 = inputdata.tofmax;
		stokes.z0 = inputdata.swl;
		stokes.gravity = inputdata.gravity;

		// set the properties of the wave
		stokes.set_stokes5_properties(stokes.wave_length, stokes.wave_height);
		if (!inputdata.property_read) {
			std::cerr << "INPUTFILE ERROR: Stokes wave selected, but no wave properties found in input file." << std::endl;
			exit(-1);
		}
		else {
			std::cout << "Stokes 5th wave initialized and ready to go." << std::endl;
		}

	}
	

	if (inputdata.wavetype == 4) {
		lsgrid.water_depth = inputdata.depth;
		lsgrid.swl = inputdata.swl;
		//lsgrid.set_ignore();

		if (lsgrid.ignore_at_init) {
			lsgrid.initialize_surface_elevation_with_ignore(irreg, lsgrid.t0);
			lsgrid.initialize_kinematics_with_ignore(irreg);
		}
		else {
			lsgrid.initialize_surface_elevation(irreg, lsgrid.t0);
			lsgrid.initialize_kinematics(irreg);
		}
	}
	if (inputdata.wavetype == 5) {
		lsgrids.water_depth = inputdata.depth;
		lsgrids.swl = inputdata.swl;
		lsgrids.set_ignore();

		if (lsgrids.ignore_at_init) {
			lsgrids.initialize_kinematics_with_ignore(irreg);
		}
		else {
			lsgrids.initialize_kinematics(irreg);
		}
	}

	
	// swd
	if (inputdata.wavetype >= 30 && inputdata.wavetype < 40) {
#if defined(SWD_enable)
		double x0_, y0_, t0_, beta_;

		std::cout << inputdata.swdFileName.c_str() << std::endl;
		// Initialize spectral wave data
		swd = new SpectralWaveData(inputdata.swdFileName.c_str(), inputdata.x_pos, inputdata.y_pos, inputdata.tofmax, inputdata.mtheta, inputdata.rho, inputdata.nsumx, inputdata.nsumy, inputdata.impl, inputdata.ipol, inputdata.norder, inputdata.dc_bias);

		std::string cid = swd->GetChr("cid");

		std::vector<std::string> v;

		std::stringstream ss(cid);

		while (ss.good()) {
			std::string substr;
			getline(ss, substr, ',');
			v.push_back(substr);
		}


		std::cout << "************************************************************************* " << std::endl;
		std::cout << "The following info is stored in the specified swd file " << inputdata.swdFileName << std::endl;
		std::cout << "************************************************************************* " << std::endl;
		for (size_t i = 0; i < v.size(); i++)
			std::cout << v[i] << std::endl;
		std::cout << "************************************************************************* " << std::endl;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;

#endif
	}

	// swd with lsgrid
	if (inputdata.wavetype == 34) {
#if defined(SWD_enable)
		double depth_swd_file = swd->GetReal("d");
		if (depth_swd_file != inputdata.depth && depth_swd_file > 0) {
			std::cout << "Warning: Specified water depth, " << inputdata.depth << "m,  not the same as used in swd file d = " << depth_swd_file << "m. Depth specified in SWD file will be used." << std::endl;
			lsgrid.water_depth = depth_swd_file;
		}
		else {
			lsgrid.water_depth = inputdata.depth;
		}

		lsgrid.swl = inputdata.swl;
		lsgrid.set_ignore();
		if (lsgrid.ignore_at_init) {
			lsgrid.initialize_surface_elevation_with_ignore(swd, lsgrid.t0);
			lsgrid.initialize_kinematics_with_ignore(swd);
		}
		else {
			lsgrid.initialize_surface_elevation(swd, lsgrid.t0);
			lsgrid.initialize_kinematics(swd);
		}
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
#if defined(VTK_enable)
	// VTK kinematics input
	if (inputdata.wavetype == 41) {
		vtkreader.init(0.);
		inputdata.depth = -vtkreader.zmin;
	}
#endif
	std::cout << "WaveID: " << inputdata.wavetype << std::endl;

	return 0;
}


double wave_water_depth() {
	return inputdata.depth;
}

/*
----------------------------------------------------------------------
Wave kinematics functions
----------------------------------------------------------------------
*/

double wave_VeloX(double xpt, double ypt, double zpt, double tpt)
{

	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-inputdata.depth, zpt);

	switch (inputdata.wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.u(tpt, xpt, ypt, zpt);
		// irregular LSgrid waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait
			sgrid.update(irregular, tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * sgrid.u(tpt, xpt, ypt, zpt);
	case 5:
	{
		if (!sgrids.CheckTime(tpt)) {
#pragma omp single nowait
			sgrids.update(irregular, tpt);
		}
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		//std::for_each(v.begin(), v.end(), [](double &el){el *= ramp.ramp(tpt, xpt, ypt); });
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[1]*ramp.ramp(tpt, xpt, ypt);
	}

		// stokes 5th
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.u(tpt, xpt, ypt, zpt);

		// wavemaker
	case 11:
		return ramp.ramp(tpt, xpt, ypt) * wavemaker.u_piston(tpt);

		// swd
	case 31:
	{
#if defined(SWD_enable)
		// Tell the swd object current application time...
		try {
			swd->UpdateTime(tpt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
		vector_swd U = swd->GradPhi(xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * U.x;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif

	}
	case 34:
	{
#if defined(SWD_enable)
		if (!sgrid.CheckTime(tpt)) {
			//#pragma omp single
#pragma omp single nowait		
			sgrid.update(swd, tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * sgrid.u(tpt, xpt, ypt, zpt);
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt)) {
#pragma omp single
			vtkreader.update(tpt);
		}

		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.u(tpt, xpt, ypt, zpt);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	
	
	default:
		return 0.0;
	}


}

double* wave_Kinematics(double xpt, double ypt, double zpt, double tpt) {
	
	double* temp;
	
	switch (inputdata.wavetype) {
	// VTKinput
	case 5:
	{
		if (!sgrids.CheckTime(tpt)) {
#pragma omp single nowait
			sgrids.update(irregular, tpt);
		}
		std::vector<double> v;
		double* res = new double[4];
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		res[0] = v[0]*ramp.ramp(tpt, xpt, ypt);
		res[1] = v[1]*ramp.ramp(tpt, xpt, ypt);
		res[2] = v[2]*ramp.ramp(tpt, xpt, ypt);
		res[3] = v[3]*ramp.ramp(tpt, xpt, ypt);
		//double rr = ramp.ramp(tpt, xpt, ypt);
		//std::for_each(v.begin(), v.end(), [](double &el){el *= rr; });
		//std::transform(v.begin(), v.end(), v.begin(),[&rr](double element) { return element *= rr; });
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		//double* a = &v[0];
		
		return res;
	}
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt)) {
#pragma omp single
			vtkreader.update(tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		double res[5];
		double* temp = vtkreader.trilinear_interpolation(res, tpt, xpt, ypt, zpt);
		return temp;
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}


	default:
		return temp;
	}
}


double wave_VeloY(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-inputdata.depth, zpt);

	switch (inputdata.wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.v(tpt, xpt, ypt, zpt);
		// irregular gridded waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.v(tpt, xpt, ypt, zpt);
	case 5:
	{
		if (!sgrids.CheckTime(tpt)) {
#pragma omp single nowait
			sgrids.update(irregular, tpt);
		}
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[2]*ramp.ramp(tpt, xpt, ypt);
	}
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.v(tpt, xpt, ypt, zpt);
		// swd
	case 31:
	{
#if defined(SWD_enable)
		// Tell the swd object current application time...
		try {
			swd->UpdateTime(tpt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
		vector_swd U = swd->GradPhi(xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * U.y;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// swd + lsgrid
	case 34:
	{
#if defined(SWD_enable)
		if (!sgrid.CheckTime(tpt)) {
			//#pragma omp single
#pragma omp single nowait
			sgrid.update(swd, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.v(tpt, xpt, ypt, zpt);
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt)) {
#pragma omp single
			vtkreader.update(tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.v(tpt, xpt, ypt, zpt);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		return 0.0;
	}
}


double wave_VeloZ(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-inputdata.depth, zpt);

	switch (inputdata.wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.w(tpt, xpt, ypt, zpt);
		// irregular gridded waves
	case 4:
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.w(tpt, xpt, ypt, zpt);
	case 5:
	{
		if (!sgrids.CheckTime(tpt)) {
#pragma omp single nowait
			sgrids.update(irregular, tpt);
		}
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[3]*ramp.ramp(tpt, xpt, ypt);
	}
	case 21:
		return ramp.ramp(tpt, xpt, ypt) * stokes5.w(tpt, xpt, ypt, zpt);
		// swd
	case 31:
	{
#if defined(SWD_enable)
		// Tell the swd object current application time...
		try {
			swd->UpdateTime(tpt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
		vector_swd U = swd->GradPhi(xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * U.z;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// Swd + lsgrid
	case 34:
	{
#if defined(SWD_enable)
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single
			sgrid.update(swd, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.w(tpt, xpt, ypt, zpt);
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt)) {
#pragma omp single
			vtkreader.update(tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.w(tpt, xpt, ypt, zpt);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		return 0.0;
	}
}

// todo: not fully functional, but not really used either.... should be updated at some point
double wave_DynPres(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	zpt = std::max(-inputdata.depth, zpt);

	switch (inputdata.wavetype) {
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
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
	{
		update_probes(tpt);
		return ramp.ramp(tpt, xpt, ypt) * irregular.eta(tpt, xpt, ypt);
	}
	case 4: {
		update_probes(tpt);
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait
			sgrid.update(irregular, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.eta(tpt, xpt, ypt);}
	case 5:
	{
		update_probes(tpt);
		if (!sgrids.CheckTime(tpt)) {
#pragma omp single nowait
			sgrids.update(irregular, tpt);
		}
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, 0., inputdata.depth);
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[0]*ramp.ramp(tpt, xpt, ypt);
	}
	case 11:
	{
		update_probes(tpt);
		return ramp.ramp(tpt, xpt, ypt) * wavemaker.wave_elev_piston(tpt);
	}
	case 21:
	{
		update_probes(tpt);
		return ramp.ramp(tpt, xpt, ypt) * stokes5.eta(tpt, xpt, ypt);
	}
	case 31:
	{
#if defined(SWD_enable)
		// Tell the swd object current application time...
		try {
			swd->UpdateTime(tpt);
		}
		catch (SwdInputValueException& e) {  //Could be t > tmax from file.
			std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
			// If we will try again with a new value of t
			// we first need to call: swd.ExceptionClear()
			exit(EXIT_FAILURE);  // In this case we just abort.
		}
		//std::cout << "time: " << tpt << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * swd->Elev(xpt, ypt);
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
		exit(-1);
#endif
	}
	case 34:
	{
#if defined(SWD_enable)
		if (!sgrid.CheckTime(tpt)) {
#pragma omp single nowait 	
			sgrid.update(swd, tpt);
		}
		return ramp.ramp(tpt, xpt, ypt) * sgrid.eta(tpt, xpt, ypt);
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt)) {
#pragma omp single
			vtkreader.update(tpt);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.eta(tpt, xpt, ypt);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		return 0.0;
	}
}

double wave_Seabed(double xpt, double ypt)
{
	switch (inputdata.wavetype) {
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		return vtkreader.seabed(xpt, ypt);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		return -inputdata.depth;;
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




void update_probes(double tpt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1: {
		if (probes.checkTime(tpt)) {
#pragma omp single nowait
			probes.write(tpt, irregular, ramp);
		}
		break;
	}
	case 4: {
#pragma omp single nowait
		if (probes.checkTime(tpt)) {
			probes.write(tpt, sgrid, irregular, ramp);
		}
		break;
	}
	case 5: {
#pragma omp single nowait
		if (probes.checkTime(tpt)) {
			probes.write(tpt, sgrids, irregular, ramp);
		}
		break;
	}
	case 11: {
		if (probes.checkTime(tpt)) {
#pragma omp single nowait
			probes.write(tpt, wavemaker, ramp);
		}
		break;
	}
	case 21: {
		if (probes.checkTime(tpt)) {
#pragma omp single nowait
			probes.write(tpt, stokes5, ramp);
		}
		break;
	}
#if defined(SWD_enable)
	case 31:
	{
		if (!probes.checkTime(tpt)) {
#pragma omp single nowait
			probes.write(tpt, swd, ramp);
		}
		break;
	}
	case 34:
	{
		if (!probes.checkTime(tpt)) {
#pragma omp single nowait
			probes.write(tpt, sgrid, swd, ramp);
		}
		break;
	}
#endif	
	}
}



//EXPORT int Init(double& tmin_in, double& tmax_in)
int wave_Initialize()
{
	std::string lineA;
	std::ifstream fid;
	std::string res;
	std::cout << "\n\n***********************************************\n\n" << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	std::cout << "CFD WAVEMAKER v3.0.1" << std::endl;
	std::cout << "---------------------------------------" << std::endl;
	
	
	//std::filesystem::path cwd = std::filesystem::current_path();

	std::string cwd = get_current_dir();

	std::cout <<"working directory: " << cwd << std::endl;

	// Check for waveinput.dat file in most common locations

	// Open and read file to find which input file version to read.
	// READ INPUT FILE AND REMOVE COMMENT LINES
	fid.open("./waveinput.dat");
	// check one step up in the folder tree (this is used in the latest comflow version)
	if (fid.fail()) {
		fid.open("./constant/waveinput.dat");
		if (!fid.fail()) {
			std::cout << "Waveinput.dat file found in folder ./constant/" << std::endl;
		}
	}
	// Special cases for comflow
	if (fid.fail()) {
		fid.open("../input_files/waveinput.dat");
		if (!fid.fail()) {
			std::cout << "Waveinput.dat file found in folder ./input_files/" << std::endl;
		}
	}
	if (fid.fail()) {
		fid.open("../waveinput.dat");
		if (!fid.fail()) {
			std::cout << "Waveinput.dat file found in folder ../" << std::endl;
		}
	}
	if (fid.fail()) {
		std::cerr << "Could not open file (is it really there?) " << std::endl;
		return -1;
		exit(1);
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
	
	//check file version
	//while (!f.eof()) {
	// Read first line only
	getline(f, lineA);
	trim(lineA);
	std::cout << lineA.substr(2, 3) << std::endl;
	std::cout << lineA.substr(0, 2).compare("@v") << std::endl;
	if (lineA.substr(0,2).compare("@v")) {	
		
		std::cout << "Please specified version at the first line of input file (hint: @v213)" << std::endl;
		exit(0);
	}
	if (stoi(lineA.substr(2,3)) < 213 ) {
		std::cout << "Your input file is too old. Please update according to the latest manual of CFDwavemaker." << std::endl;
		exit(0);
	}


	//}

	

	int i = process_inputdata(res, irregular, stokes5, wavemaker, sgrid, sgrids, ramp);
	CFDwmInit = true;
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

int CFDwavemaker_is_initialized() {
	return int(CFDwmInit);
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
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return irregular.phase_velocity(opt);
	case 4:
		return irregular.phase_velocity(opt);
	default:
		return 0.;
	}
}

double wave_mean_length(int opt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return irregular.mean_wave_length(opt);
	case 4:
		return irregular.mean_wave_length(opt);
	default:
		return 0.;
	}
}
double wave_mean_period(int opt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return irregular.mean_wave_period(opt);
	case 4:
		return irregular.mean_wave_period(opt);
	default:
		return 0.;
	}	
}

void wave_force_update(double tpt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 4: {
		if (!sgrid.CheckTime(tpt)) {
			//#pragma omp single
#pragma omp single nowait
			sgrid.update(irregular, tpt);
			}
			break;
		}
	case 5: {
		if (!sgrids.CheckTime(tpt)) {
			//#pragma omp single
#pragma omp single nowait
			sgrids.update(irregular, tpt);
		}
		break;
	}
#if defined(SWD_enable)
	case 34: {
		
		if (!sgrid.CheckTime(tpt)) {
			//#pragma omp single
#pragma omp single nowait 			
			sgrid.update(swd, tpt);
		}
		break;
	}
#endif

#if defined(VTK_enable)
	case 41: {
		if (!vtkreader.CheckTime(tpt)) {
#pragma omp single
			vtkreader.update(tpt);
		}
		break;
	}
#endif
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