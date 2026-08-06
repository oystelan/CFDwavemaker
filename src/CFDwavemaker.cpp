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
// Version: set in src/VERSION (passed in by the makefiles as CFDWM_VERSION)
#ifndef CFDWM_VERSION
#define CFDWM_VERSION "(unversioned build)"
#endif
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
#include <math.h>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort
//#include <filesystem>


//#include <direct.h> // windows only function
#//include <cctype>
//#include <locale>
#include "CFDwavemaker.h"
#include "Stokes5.h"
#include "FentonStream.h"
#include "Irregular.h"
#include "Utils.h"
#include "Wavemaker.h"
#include "lsgrid_spline.h"
#include "mpi_utils.h"
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
	double x_pos, y_pos, tofmax;
	double mtheta = 0.;
	double swl = 0.;
	bool bw_auto_calc = false;
	int wavetype;
	bool property_read = false;
	double rho = 1.;
	double gravity = 9.81;
	double ampl = 1.;
	// SWD parameters
	double nsumx = -1, nsumy = -1, impl = 0, ipol = 0, norder = 0;
	bool dc_bias = false;
	std::string swdFileName;

	// sloping water related parameters
	double slope_angle = 0.; // tilt angle of the slope
	double slope_direction = 0.; // direction of water slope (degrees)
	double slope_rotmat[9] = {}; // store rotation matrix (flattened)
	double slope_radius_limiter = HUGE;
	// lsgrid related data
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

	// Nonhydrostatic pressure options
	// When false (default): p_nh = -rho*(dphi/dt + g*eta)  [linearized, suitable for non-hydrostatic CFD solvers]
	// When true:            p_nh = -rho*(dphi/dt + 0.5*|u|^2 + g*eta)  [full Bernoulli form]
	bool nhpres_include_kinetic = false;

	// Regular wave (wavetype 21) theory selection under [wave input regular]:
	//   0 = Stokes 5th (default), 1 = Fenton stream function.
	int  regular_theory = 0;
	int  wave_order = 16;                 // Fourier order N for stream function
	bool regular_specify_length = true;   // true: user gave wave_length; false: wave_period

#if defined(SWD_enable)
	// Spectral properties derived from the SWD file at init (long-crested shapes only)
	SwdSpectralProps swd_spectral;
#endif

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

// Fenton stream function class (regular waves, theory_type stream)
FentonStream fenton;

// Irregular class
Irregular irregular;

// Wavemaker theory class
Wavemaker wavemaker;

// Grid class
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
int process_inputdata(std::string res, Irregular& irreg, Stokes5& stokes, Wavemaker& wmaker, lsGridSpline& lsgrids, Ramp& rramp) {
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
		// others: 51-60
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
			else if (!lineA.compare("piston")) {
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
			else if (!lineA.compare("slope")) {
				inputdata.wavetype = 51;
				std::cout << "Slope of water. Water at rest" << std::endl;
			}
			
			else {
				std::cerr << "INPUTFILE ERROR: Unknown wave type specified. Valid alternatives are:" << std::endl;
				std::cerr << "irregular" << std::endl;
				std::cerr << "regular" << std::endl;
				std::cerr << "piston" << std::endl;
				std::cerr << "swd" << std::endl;
				std::cerr << "vtk" << std::endl;
				std::cerr << "slope" << std::endl;
				
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
					std::cout << "Mean wave direction: " << inputdata.mtheta << " degrees" << std::endl;
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
					buf >> inputdata.ampl;
					buf.clear();
					std::cout << "Amplify (gain): " << inputdata.ampl << "m" << std::endl;
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
					std::cout << "Water density rho set to: " << inputdata.rho << "kg/m^3" << std::endl;
				}
				if (!lineA.compare(0, 22, "nhpres_include_kinetic")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;
					if (!dummystr.compare(0, 4, "true")) {
						inputdata.nhpres_include_kinetic = true;
						std::cout << "Nonhydrostatic pressure: full Bernoulli form (incl. kinetic energy term)" << std::endl;
					}
					else {
						inputdata.nhpres_include_kinetic = false;
						std::cout << "Nonhydrostatic pressure: linearized form (kinetic energy term excluded)" << std::endl;
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
					rramp.ramp_init_time_up = true;
					buf >> rramp.time_rampup_start;
					buf >> rramp.time_rampup_end;
					buf.clear();
					std::cout << "time-rampup -- starttime: " << rramp.time_rampup_start  << " sec, endtime: " << rramp.time_rampup_end << " sec." << std::endl;
				}
				if (!lineA.compare(0, 13, "time_rampdown")) {
					buf.str(lineA);
					buf >> dummystr;
					rramp.ramp_init_time_down = true;
					buf >> rramp.time_rampdown_start;
					buf >> rramp.time_rampdown_end;
					buf.clear();
					std::cout << "time-rampdown -- starttime: " << rramp.time_rampdown_start << " sec, endtime: " << rramp.time_rampdown_end << " sec." << std::endl;
				}
				if (!lineA.compare(0, 8, "x_rampup")) {
					buf.str(lineA);
					buf >> dummystr;
					rramp.ramp_init_x_up = true;
					buf >> rramp.x_rampup_start;
					buf >> rramp.x_rampup_end;
					buf.clear();
					std::cout << "x-rampup -- startpos: " << rramp.x_rampup_start << " m, endpos: " << rramp.x_rampup_end << " m." << std::endl;
				}
				if (!lineA.compare(0, 10, "x_rampdown")) {
					buf.str(lineA);
					buf >> dummystr;
					rramp.ramp_init_x_down = true;
					buf >> rramp.x_rampdown_start;
					buf >> rramp.x_rampdown_end;
					buf.clear();
					std::cout << "x-rampdown -- startpos: " << rramp.x_rampdown_start << " m, endpos: " << rramp.x_rampdown_end << " m." << std::endl;
				}
				if (!lineA.compare(0, 8, "y_rampup")) {
					buf.str(lineA);
					buf >> dummystr;
					rramp.ramp_init_y_up = true;
					buf >> rramp.y_rampup_start;
					buf >> rramp.y_rampup_end;
					buf.clear();
					std::cout << "y-rampup -- startpos: " << rramp.y_rampup_start << " m, endpos: " << rramp.y_rampup_end << " m." << std::endl;
				}
				if (!lineA.compare(0, 10, "y_rampdown")) {
					buf.str(lineA);
					buf >> dummystr;
					rramp.ramp_init_y_down = true;
					buf >> rramp.y_rampdown_start;
					buf >> rramp.y_rampdown_end;
					buf.clear();
					std::cout << "y-rampdown -- startpos: " << rramp.y_rampdown_start << " m, endpos: " << rramp.y_rampdown_end << " m." << std::endl;
				}
				
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}
		// Wave properties: irregular wave components (manual)
		if (!lineA.compare("[wave input irregular]")) {
			std::cout << "--------------------------" << std::endl;
			std::cout << "Irregular input data:" << std::endl;
			std::cout << "--------------------------" << std::endl;

			if (inputdata.wavetype != 1) {
				std::cerr << "INPUTFILE ERROR: irregular wave components does not match the specified wave type. Check inputfile" << std::endl;
				exit(-1);
			}
			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				//std::cout << lineA << std::endl;
				if (!lineA.compare(0, 22, "kinematics_surf_extrap")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> irreg.extrapolation_met;
					buf.clear();
					if (irreg.extrapolation_met == 2){
						irreg.order = 2;
						std::cout << "Second order kinematics used with Taylor expansion above SWL. " << std::endl;
					}
					else if (irreg.extrapolation_met == 1){
						irreg.order = 0;
						std::cout << "Linear wave theory used with constant value above swl." << std::endl;
					}
					else if (irreg.extrapolation_met == 0){
						irreg.order = 0;
						std::cout << "Linear wave theory used with exponential extrapolation above z=swl." << std::endl;
					}
					else{
						std::cout << "Unknown extrapolation method. try a value between 0 and 2." << std::endl;
						exit(-1);
					}

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
				if (!lineA.compare(0, 15, "individual_ramp")) {

					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;
					if (!dummystr.compare(0, 4, "true")) {
						std::cout << "Individual rampup of each frequency component set to true. Three additional columns of data required when specifying frequency components." << std::endl;
						irreg.individual_ramp = true;
					}
					buf.clear();
				}
				if (!lineA.compare(0, 15, "vertical_theory")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;
					if (!dummystr.compare(0, 4, "smit")) {
						irreg.vertical_theory = 1;
						std::cout << "Smit et al. (2017) s-coordinate vertical theory selected." << std::endl;
						std::cerr << "WARNING: vertical_theory 'smit' is EXPERIMENTAL / in development and not yet "
							"validated for production use (the off-diagonal second-order terms are still being "
							"verified). Use the default 'sharma' theory unless you are specifically testing this." << std::endl;
					}
					else {
						irreg.vertical_theory = 0;
						std::cout << "Conventional Sharma & Dean vertical theory (default)." << std::endl;
					}
					buf.clear();
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 5, "nfreq")) {
					break;
				}
			}
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
			if (irreg.ndir == 0 && irreg.individual_ramp) {
				std::cout << "Irregular seas, one directional component for each frequency specified." << std::endl;
				std::cout << "Rampup start time for each frequency components specified sixed column." << std::endl;
				std::cout << "Number of frequency components: " << irreg.nfreq << std::endl;
				irreg.ndir = 1;

				// Create some temporary vectors for storage of spectral data
				std::vector<double> omega;
				std::vector<double> A;
				std::vector<double> k;
				std::vector<double> theta;
				std::vector<double> phase;
				std::vector<double> freq_ramptime;
				std::vector<double> freq_starttime;
				std::vector<double> freq_duration;

				double temp;
				std::cout << "# OMEGA[rad / s]    A[m]         K[1/m],            Phase[rad],        Theta [rad],   Rampup[sec], ramp start time[sec], freq. duration[sec]" << std::endl;
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
					buf >> temp;
					freq_ramptime.push_back(temp);
					buf >> temp;
					freq_starttime.push_back(temp);
					buf >> temp;
					freq_duration.push_back(temp);
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
					irreg.freq_ramptime.push_back(freq_ramptime[i]);
					irreg.freq_starttime.push_back(freq_starttime[i]);
					irreg.freq_duration.push_back(freq_duration[i]);
				}
				irreg.initialized = true;

			}
			else if (irreg.ndir == 0) {
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
					irreg.freq_ramptime.push_back(0.);
					irreg.freq_duration.push_back(HUGE);
					irreg.freq_starttime.push_back(-SMALL);
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
					irreg.freq_ramptime.push_back(0.);
					irreg.freq_duration.push_back(HUGE);
					irreg.freq_starttime.push_back(-SMALL);
				}
				irreg.initialized = true;
			}
			inputdata.property_read = true;
		}

		// Wave properties: piston wave maker
		if (!lineA.compare("[wave input piston]")) {

			if (inputdata.wavetype != 11) {
				std::cerr << "piston wave maker properties does not match the specified wave type. Check inputfile" << std::endl;
				exit(-1);
			}
			if (wmaker.initialized) {
				std::cerr << "INPUTFILE ERROR: Pistonwavemaker already initialized. please check input file" << std::endl;
				exit(-1);
			}
			// Wavemaker theory (piston)
			getline(f, lineA);
			buf.str(lineA);
			buf >> wmaker.n_timesteps;
			//n_timesteps = stoi(lineA);
			std::cout << "Number of timesteps: " << wmaker.n_timesteps << std::endl;

			// declare some vectors to store piston data
			wmaker.PD_time = new double[wmaker.n_timesteps];
			wmaker.PD_ampl = new double[wmaker.n_timesteps];
			wmaker.PD_velo = new double[wmaker.n_timesteps];
			wmaker.PD_eta = new double[wmaker.n_timesteps];
			wmaker.PD_velo_grad = new double[wmaker.n_timesteps];

			for (int i = 0; i < wmaker.n_timesteps; i++) {
				getline(f, lineA);
				buf.str(lineA);
				buf >> wmaker.PD_time[i];
				buf >> wmaker.PD_ampl[i];
				buf >> wmaker.PD_velo[i];
				buf >> wmaker.PD_eta[i];
				buf.clear();
			}
			wmaker.gradient1D(wmaker.PD_time, wmaker.PD_velo, wmaker.PD_velo_grad);
			inputdata.property_read = true;
			wmaker.initialized = true;
			std::cout << "piston timeseries read successfully." << std::endl;
		}
		
		if (!lineA.compare("[sloping water]")) {
			while (!f.eof()) {
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 15, "slope_direction")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.slope_direction;
					buf.clear();
					std::cout << "Direction of slope (degrees):  " << inputdata.slope_direction << std::endl;
				}
				if (!lineA.compare(0, 11, "slope_angle")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.slope_angle;
					buf.clear();
					std::cout << "Slope angle (degrees) positive value gives upward slope:  " << inputdata.slope_angle << std::endl;
				}
				if (!lineA.compare(0, 15, "limiting_radius")) {
					//  - Dict input data									
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.slope_radius_limiter;
					buf.clear();
					std::cout << "Limiting radius of slope (flat outside this radius):  " << inputdata.slope_radius_limiter << std::endl;
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}	
		}

		// Wave properties: Stokes wave properties
		if (!lineA.compare("[wave input regular]")) {
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
				if (!lineA.compare(0, 11, "theory_type")) {
					buf.str(lineA);
					buf >> dummystr;
					std::string tt;
					buf >> tt;
					buf.clear();
					if (!tt.compare("stream") || !tt.compare("fenton")) {
						inputdata.regular_theory = 1;
						std::cout << "Regular wave theory: Fenton stream function" << std::endl;
					}
					else {
						inputdata.regular_theory = 0;
						std::cout << "Regular wave theory: Stokes 5th order" << std::endl;
					}
				}
				if (!lineA.compare(0, 10, "wave_order")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> inputdata.wave_order;
					buf.clear();
					std::cout << "Stream function Fourier order: " << inputdata.wave_order << std::endl;
				}
				if (!lineA.compare(0, 11, "wave_length")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> stokes.wave_length;
					buf.clear();
					inputdata.regular_specify_length = true;
					// Make sure depth/gravity are propagated so the dispersion solver works
					stokes.depth = inputdata.depth;
					stokes.gravity = inputdata.gravity;
					stokes.wave_period = stokes.wavePeriodFromLength(stokes.wave_length);
					std::cout << "Wave length:  " << stokes.wave_length << std::endl;
					std::cout << "Wave period computed from specified wave length: " << stokes.wave_period << std::endl;
				}
				if (!lineA.compare(0, 11, "wave_period")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> stokes.wave_period;
					buf.clear();
					inputdata.regular_specify_length = false;
					stokes.depth = inputdata.depth;
					stokes.gravity = inputdata.gravity;
					stokes.wave_length = stokes.wavelengthFromPeriod(stokes.wave_period);
					std::cout << "Wave period: " << stokes.wave_period << std::endl;
					std::cout << "Wave length computed from specified wave period: " << stokes.wave_length << std::endl;
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

		if (!lineA.compare("[wave input swd]")) { //optional
			std::cout << "-------------------------------------" << std::endl;
			std::cout << "Spectral wave data:" << std::endl;
			std::cout << "-------------------------------------" << std::endl;
			if (inputdata.wavetype != 31 && inputdata.wavetype != 34) {
				std::cerr << "INPUTFILE ERROR: please set [wave type] = swd when using keyword [wave input swd]." << std::endl;
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
				if (!lineA.compare(0, 7, "dc_bias")) {
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
			

			std::cout << "Interpolation Scheme: spline" << std::endl;
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

			if (inputdata.wavetype == 1){
				inputdata.wavetype = 5;
				std::cout << "LS grid + Spline + irregular" << std::endl;
			}
			else if (inputdata.wavetype == 31){
				inputdata.wavetype = 34;
				std::cout << "LS grid + Spline + swd" << std::endl;
			}
			else {
				std::cerr << "INPUTFILE ERROR: LS grid may currently only be used with irregular second order wave theory or swd" << std::endl;
				exit(-1);
			}
			// copy necessary parameters to lsgrid spline class
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
		
		if (!lineA.compare("[lsgrid2vtk]")) { //optional
			std::cout << "--------------------" << std::endl;
			std::cout << "VTK output settings:" << std::endl;
			std::cout << "--------------------" << std::endl;
			if (inputdata.wavetype != 5)  {
				std::cout << "InputError: VTK output only works in combination with LSgrid at the moment." << std::endl;
				exit(-1);
			}
			sgrids.dump_vtk = true;

			while (!f.eof()) {
				lineP = lineA;
				getline(f, lineA);
				trim(lineA);
				if (!lineA.compare(0, 12, "storage_path")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> lsgrids.vtk_directory_path;
					buf.clear();
					std::cout << "Directory for storage of vtk files: " << lsgrids.vtk_directory_path << std::endl;
				}
				if (!lineA.compare(0, 8, "filename")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> lsgrids.vtk_prefix;
					buf.clear();
					std::cout << "filename prefix: " << lsgrids.vtk_prefix << std::endl;
				}
				if (!lineA.compare(0, 13, "vtk_timelabel")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> lsgrids.vtk_timelabel;
					buf.clear();
					std::cout << "time label to use in vtk specified as: " << lsgrids.vtk_timelabel << std::endl;
				}
				// if new tag is reach. break while loop.
				if (!lineA.compare(0, 1, "[")) {
					skip_getline = true;
					break;
				}
			}
		}


		if (!lineA.compare("[wave input vtk]")) {

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
				// Multilayer vertical velocity reconstruction (only used when cell field 'h' is present):
				//   "ppr" (default, conservative parabolic) or "linear" (piecewise linear between cell centres).
				if (!lineA.compare(0, 15, "vertical_interp")) {
					buf.str(lineA);
					buf >> dummystr;
					buf >> dummystr;
					vtkreader.vertical_interp = (!dummystr.compare(0, 6, "linear")) ? 1 : 0;
					std::cout << "Multilayer vertical velocity reconstruction: "
						<< (vtkreader.vertical_interp == 0 ? "PPR" : "linear") << std::endl;
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

	if (inputdata.wavetype == 5 && irregular.order != 2) {
		std::cerr << "LSgrid interpolation uses strictly second order theory. Doesnt make sense to use grid interpolation for linear theory since this is very fast anyway. Set kinematics_surf_extrap to 2 under [wave input irregular], or remove [lsgrid]." << std::endl;
		exit(-1);
	}


	std::cout << "\n-----------------------------------------------" << std::endl;
	std::cout << "Input file read successfully." << std::endl;
	std::cout << "***********************************************\n\n" << std::endl;

	std::cout << "Wavetype: " << inputdata.wavetype << std::endl;


	if (inputdata.wavetype == 1 || inputdata.wavetype == 5) {
		irreg.depth = inputdata.depth;
		irreg.mtheta = inputdata.mtheta;
		irreg.tofmax = inputdata.tofmax;
		irreg.fpoint[0] = inputdata.x_pos;
		irreg.fpoint[1] = inputdata.y_pos;
		irreg.swl = inputdata.swl;
		irreg.g = inputdata.gravity;
		irreg.rho = inputdata.rho;
		irreg.ampl = inputdata.ampl;

		if (inputdata.bw_auto_calc) {
			irreg.dw_bandwidth = irregular.bandwidth_estimator();
			std::cout << "BW_autocalc: Bandwidth parameter set to " << irreg.dw_bandwidth << " rad/s." << std::endl;		
		}
		irreg.normalize_data();
		irreg.calculate_bwindices();
		irreg.dumpSpectralComponents();
		irreg.print_spectral_parameters();
		std::cout << "irregorder: " << irreg.order << std::endl;
		if (irreg.order == 2){
			std::cout << "done initializing" << std::endl;
			irreg.rebuild_second_order_cache();
			std::cout << "done initializing" << std::endl;
		}

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
		if (!inputdata.property_read) {
			std::cerr << "INPUTFILE ERROR: Regular wave selected, but no wave properties found in input file." << std::endl;
			exit(-1);
		}
		if (inputdata.regular_theory == 1) {
			// Fenton stream function wave
			fenton.depth = inputdata.depth;
			fenton.theta = inputdata.mtheta;
			fenton.x0 = inputdata.x_pos;
			fenton.y0 = inputdata.y_pos;
			fenton.t0 = inputdata.tofmax;
			fenton.z0 = inputdata.swl;
			fenton.gravity = inputdata.gravity;
			fenton.rho = inputdata.rho;
			fenton.current = stokes.current;   // shares the current_speed keyword

			// Solve: specify_length uses wave_length, else wave_period. Current is Eulerian mean.
			fenton.set_stream_properties(stokes.wave_height, stokes.wave_length, stokes.wave_period,
				inputdata.regular_specify_length, /*current_type=Eulerian*/0,
				inputdata.wave_order, /*height-continuation steps*/5);
			std::cout << "Fenton stream function wave initialized (N=" << inputdata.wave_order
				<< ", L=" << fenton.wave_length << ", T=" << fenton.wave_period
				<< ", c=" << fenton.c << ")." << std::endl;
		}
		else {
			// Stokes 5th order
			stokes.depth = inputdata.depth;
			stokes.theta = inputdata.mtheta;
			stokes.x0 = inputdata.x_pos;
			stokes.y0 = inputdata.y_pos;
			stokes.t0 = inputdata.tofmax;
			stokes.z0 = inputdata.swl;
			stokes.gravity = inputdata.gravity;
			stokes.rho = inputdata.rho;

			stokes.set_stokes5_properties(stokes.wave_length, stokes.wave_height);
			std::cout << "Stokes 5th wave initialized and ready to go." << std::endl;
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

		// Compute spectral wave properties from the SWD file (long-crested shapes only)
		inputdata.swd_spectral = compute_swd_spectral_properties(swd, inputdata.gravity, inputdata.depth);
		if (inputdata.swd_spectral.valid) {
			std::cout << "Spectral properties derived from SWD file:" << std::endl;
			std::cout << "  Hs = " << inputdata.swd_spectral.Hs << " m" << std::endl;
			std::cout << "  Tp = " << inputdata.swd_spectral.Tp << " s (peak period)" << std::endl;
			std::cout << "  T1 = " << inputdata.swd_spectral.T1 << " s (first-moment period, m0/m1)" << std::endl;
			std::cout << "  T2 = " << inputdata.swd_spectral.T2 << " s (zero-crossing period, sqrt(m0/m2))" << std::endl;
			std::cout << "************************************************************************* " << std::endl;
		}
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;

#endif
	}

	// SWD + lsgrid (case 34): sample the swd kinematics onto the spline grid.
	// Must come after the swd object is created above. The fill is serial
	// within the process (the swd library allows only one object per process
	// and is stateful); MPI ranks each fill their own column block.
	if (inputdata.wavetype == 34) {
#if defined(SWD_enable)
		lsgrids.water_depth = inputdata.depth;
		lsgrids.swl = inputdata.swl;
		lsgrids.set_ignore();

		if (lsgrids.ignore_at_init) {
			lsgrids.initialize_kinematics_with_ignore(swd);
		}
		else {
			lsgrids.initialize_kinematics(swd);
		}
#endif
	}

#if defined(VTK_enable)
	// VTK kinematics input
	if (inputdata.wavetype == 41) {
		vtkreader.init(0.);
		inputdata.depth = -vtkreader.zmin;
	}
#endif
	if (inputdata.wavetype == 51 ) {

		double yaw = inputdata.slope_direction*PI/180.;
		double pitch = inputdata.slope_angle*PI/180.;
		// Sloping water at rest. Compute rotation matrix from given angles for slope
		inputdata.slope_rotmat[0] = cos(yaw)*cos(pitch)*cos(yaw)+sin(yaw)*sin(yaw);
		inputdata.slope_rotmat[1] = cos(yaw)*cos(pitch)*sin(yaw)-sin(yaw)*cos(yaw);
		inputdata.slope_rotmat[2] = -cos(yaw)*sin(pitch);
		inputdata.slope_rotmat[3] = sin(yaw)*cos(pitch)*cos(yaw)-cos(yaw)*sin(yaw);
		inputdata.slope_rotmat[4] = sin(yaw)*cos(pitch)*sin(yaw)+cos(yaw)*cos(yaw);
		inputdata.slope_rotmat[5] = -sin(yaw)*sin(pitch);
		inputdata.slope_rotmat[6] = sin(pitch)*cos(yaw);
		inputdata.slope_rotmat[7] = sin(pitch)*sin(yaw);
		inputdata.slope_rotmat[8] = cos(pitch);

	}
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
	case 34: // swd + lsgrid: served from the same spline grid
	case 5:
	{
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		//std::for_each(v.begin(), v.end(), [](double &el){el *= ramp.ramp(tpt, xpt, ypt); });
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[1]*ramp.ramp(tpt, xpt, ypt);
	}

		// regular waves (Stokes 5th / Fenton stream function)
	case 21:
		if (inputdata.regular_theory == 1)
			return ramp.ramp(tpt, xpt, ypt) * fenton.u(tpt, xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * stokes5.u(tpt, xpt, ypt, zpt);

		// wavemaker
	case 11:
		return ramp.ramp(tpt, xpt, ypt) * wavemaker.u_piston(tpt);

		// swd
	case 31:
	{
#if defined(SWD_enable)
		// SWD has a stateful API (UpdateTime mutates internal state).
		// Serialize all SWD access across threads via named critical section.
		double ux;
#pragma omp critical(swd_section)
		{
			try {
				swd->UpdateTime(tpt);
			}
			catch (SwdInputValueException& e) {  //Could be t > tmax from file.
				std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}
			double eta_val = swd->Elev(xpt, ypt);
			vector_swd U = swd->GradPhi(xpt, ypt, std::min(zpt,eta_val));
			ux = U.x;
		}
		return ramp.ramp(tpt, xpt, ypt) * ux;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif

	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}

		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.u(tpt - inputdata.tofmax, xpt - inputdata.x_pos, ypt - inputdata.y_pos, zpt);
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
	case 34: // swd + lsgrid: served from the same spline grid
	case 5:
	{
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
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		double res[5];
		double* temp = vtkreader.trilinear_interpolation(res, tpt - inputdata.tofmax, xpt - inputdata.x_pos, ypt - inputdata.y_pos, zpt);
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
	case 34: // swd + lsgrid: served from the same spline grid
	case 5:
	{
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[2]*ramp.ramp(tpt, xpt, ypt);
	}
	case 21:
		if (inputdata.regular_theory == 1)
			return ramp.ramp(tpt, xpt, ypt) * fenton.v(tpt, xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * stokes5.v(tpt, xpt, ypt, zpt);
		// swd
	case 31:
	{
#if defined(SWD_enable)
		double uy;
#pragma omp critical(swd_section)
		{
			try {
				swd->UpdateTime(tpt);
			}
			catch (SwdInputValueException& e) {  //Could be t > tmax from file.
				std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}
			double eta_val = swd->Elev(xpt, ypt);
			vector_swd U = swd->GradPhi(xpt, ypt, std::min(zpt,eta_val));
			uy = U.y;
		}
		return ramp.ramp(tpt, xpt, ypt) * uy;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.v(tpt - inputdata.tofmax, xpt - inputdata.x_pos, ypt - inputdata.y_pos, zpt);
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
	case 34: // swd + lsgrid: served from the same spline grid
	case 5:
	{
		std::vector<double> v;
		v = sgrids.get_kinematics_at_point(tpt, xpt, ypt, zpt, inputdata.depth);
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return v[3]*ramp.ramp(tpt, xpt, ypt);
	}
	case 21:
		if (inputdata.regular_theory == 1)
			return ramp.ramp(tpt, xpt, ypt) * fenton.w(tpt, xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * stokes5.w(tpt, xpt, ypt, zpt);
		// swd
	case 31:
	{
#if defined(SWD_enable)
		double uz;
#pragma omp critical(swd_section)
		{
			try {
				swd->UpdateTime(tpt);
			}
			catch (SwdInputValueException& e) {  //Could be t > tmax from file.
				std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}
			double eta_val = swd->Elev(xpt, ypt);
			vector_swd U = swd->GradPhi(xpt, ypt, std::min(zpt,eta_val));
			uz = U.z;
		}
		return ramp.ramp(tpt, xpt, ypt) * uz;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}
		
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.w(tpt - inputdata.tofmax, xpt - inputdata.x_pos, ypt - inputdata.y_pos, zpt);
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

	default:
		return 0.0;
	}
}

// Nonhydrostatic pressure: p_nh = p_total - rho*g*(eta - z)
double wave_nhPres(double xpt, double ypt, double zpt, double tpt)
{
	// Clamp to seabed
	zpt = std::max(-inputdata.depth, zpt);

	switch (inputdata.wavetype) {
	// irregular waves
	// p_nh = dp - rho*g*eta, where dp = -rho*d(phi)/dt is what irregular.dp() returns.
	// Equivalent to case 31's -rho*(dphi/dt + g*eta) (same sign convention).
	// Note: at second order, irregular.dp() includes kinetic energy contributions
	// (via dp2). The nhpres_include_kinetic flag currently only affects case 31.
	case 1:
	{
		double dp_raw = irregular.dp(tpt, xpt, ypt, zpt);
		double eta_raw = irregular.eta(tpt, xpt, ypt);
		return ramp.ramp(tpt, xpt, ypt) * (dp_raw - irregular.rho * irregular.g * eta_raw);
	}
	// Regular waves (Stokes 5th / Fenton stream function): non-hydrostatic pressure
	// via steady Bernoulli in the wave frame.
	case 21:
		if (inputdata.regular_theory == 1)
			return ramp.ramp(tpt, xpt, ypt) * fenton.dp(tpt, xpt, ypt, zpt);
		return ramp.ramp(tpt, xpt, ypt) * stokes5.dp(tpt, xpt, ypt, zpt);
	// VTK multilayer: non-hydrostatic pressure phi read from the .vts cell field
	// (defined on the layer interfaces; requires the 'phi' field to be present).
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.phi(tpt - inputdata.tofmax, xpt - inputdata.x_pos, ypt - inputdata.y_pos, zpt);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
		return 0.0;
#endif
	}
	// SWD - computed from primitives: p_nh = -rho*(dphi/dt + 0.5*|u|^2 + g*eta)
	case 31:
	{
#if defined(SWD_enable)
		double p_nh = 0.0;
		bool above_surface = false;
#pragma omp critical(swd_section)
		{
			try {
				swd->UpdateTime(tpt);
			}
			catch (SwdInputValueException& e) {
				std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}
			double eta_val = swd->Elev(xpt, ypt);
			if (zpt > eta_val) {
				above_surface = true;
			} else {
				double phi_t = swd->DdtPhi(xpt, ypt, zpt);
				double kinetic = 0.0;
				if (inputdata.nhpres_include_kinetic) {
					vector_swd U = swd->GradPhi(xpt, ypt, zpt);
					kinetic = 0.5 * (U.x * U.x + U.y * U.y + U.z * U.z);
				}
				p_nh = -inputdata.rho * (phi_t + kinetic + inputdata.gravity * eta_val);
			}
		}
		if (above_surface) return 0.0;
		return ramp.ramp(tpt, xpt, ypt) * p_nh;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}
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
	case 34: // swd + lsgrid: served from the same spline grid
	case 5:
	{
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
		if (inputdata.regular_theory == 1)
			return ramp.ramp(tpt, xpt, ypt) * fenton.eta(tpt, xpt, ypt);
		return ramp.ramp(tpt, xpt, ypt) * stokes5.eta(tpt, xpt, ypt);
	}
	case 31:
	{
#if defined(SWD_enable)
		double eta_val;
#pragma omp critical(swd_section)
		{
			try {
				swd->UpdateTime(tpt);
			}
			catch (SwdInputValueException& e) {  //Could be t > tmax from file.
				std::cout << typeid(e).name() << std::endl << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}
			eta_val = swd->Elev(xpt, ypt);
		}
		return ramp.ramp(tpt, xpt, ypt) * eta_val;
#else
		std::cerr << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
		exit(-1);
#endif
	}
	// VTKinput
	case 41:
	{
#if defined(VTK_enable)
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}
		//std::cout << zpt << " u: " << sgrid.u(tpt, xpt, ypt, zpt) << std::endl;
		return ramp.ramp(tpt, xpt, ypt) * vtkreader.eta(tpt - inputdata.tofmax, xpt - inputdata.x_pos, ypt - inputdata.y_pos);
#else
		std::cerr << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	case 51:
	{
		double pointradius = sqrt(xpt*xpt+ypt*ypt);
		double xpt_new = std::min(inputdata.slope_radius_limiter, pointradius)*std::min(1.,(xpt/pointradius));
		double ypt_new = std::min(inputdata.slope_radius_limiter, pointradius)*std::min(1.,(ypt/pointradius));
		return inputdata.swl + inputdata.slope_rotmat[6]*xpt_new+ inputdata.slope_rotmat[7]*ypt_new;

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
		return vtkreader.seabed(xpt - inputdata.x_pos, ypt - inputdata.y_pos);
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
	case 34: // swd + lsgrid: served from the same spline grid
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
			if (inputdata.regular_theory == 1)
				probes.write(tpt, fenton, ramp);
			else
				probes.write(tpt, stokes5, ramp);
		}
		break;
	}
#if defined(SWD_enable)
	case 31:
	{
		if (probes.checkTime(tpt)) {
#pragma omp single nowait
			probes.write(tpt, swd, ramp);
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
	std::cout << "CFD WAVEMAKER " << CFDWM_VERSION << std::endl;
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
	if (lineA.substr(0,2).compare("@v")) {
		std::cout << "Please specify the input file version on the first line of the input file (hint: @v310)" << std::endl;
		exit(0);
	}
	if (stoi(lineA.substr(2,3)) < 300 ) {
		std::cout << "your waveinput.dat file is outdated and does not match the current version of CFDwavemaker. See manual " << CFDWM_VERSION << " or newer for the updated format" << std::endl;
		exit(0);
	}


	//}

	
	
	int i = process_inputdata(res, irregular, stokes5, wavemaker, sgrids, ramp);
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

double wave_phase_speed(int opt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return irregular.phase_speed(opt);
	default:
		return 0.;
	}
}

double wave_mean_length(int opt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return irregular.mean_wave_length(opt);
	case 21:
		return (inputdata.regular_theory == 1) ? fenton.wave_length : stokes5.wave_length;
	default:
		return 0.;
	}
}
double wave_mean_period(int opt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 1:
		return irregular.mean_wave_period(opt);
	case 21: // regular wave, return period
		return (inputdata.regular_theory == 1) ? fenton.wave_period : stokes5.wave_period;
#if defined(SWD_enable)
	case 31:
		// SWD: spectral periods derived from DFT of eta(x) at init.
		// opt: 1 = Tp (peak), 2 = T2 (zero-crossing), otherwise T1 (first-moment, default)
		if (!inputdata.swd_spectral.valid) return 0.0;
		if (opt == 1) return inputdata.swd_spectral.Tp;
		if (opt == 2) return inputdata.swd_spectral.T2;
		return inputdata.swd_spectral.T1;
#endif
	default:
		return 0.;
	}
}

// Optional: set the MPI communicator for the distributed lsGridSpline update.
// Accepts a Fortran integer communicator handle so the public header stays
// MPI-free. No-op unless built with -DMPI_enable=1.
void wave_set_mpi_comm(int fortran_comm) {
#if defined(MPI_enable)
	cfdwm_mpi::set_comm(MPI_Comm_f2c((MPI_Fint)fortran_comm));
#else
	(void)fortran_comm;
#endif
}

void wave_force_update(double tpt) {
	switch (inputdata.wavetype) {
	case 5: {
		// while (not if): catch up if tpt has advanced more than one grid dt
		// since the last update (each update() advances the grid one dt).
		while (sgrids.CheckTime(tpt)) {
			sgrids.update(irregular, tpt);   // all threads enter here
		}
		break;
	}
#if defined(SWD_enable)
	case 34: { // swd + lsgrid: advance the grid by resampling the swd field
		while (sgrids.CheckTime(tpt)) {
			sgrids.update(swd, tpt);
		}
		break;
	}
#endif
#if defined(VTK_enable)
	case 41: {
		if (!vtkreader.CheckTime(tpt - inputdata.tofmax)) {
#pragma omp single
			vtkreader.update(tpt - inputdata.tofmax);
		}
		break;
	}
#endif
	}
}


double wave_piston_position(double tpt) {
	switch (inputdata.wavetype) {
		// return position value of piston wavemaker 
	case 11:
		return wavemaker.position_piston(tpt);
	default:
		return 0.;
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