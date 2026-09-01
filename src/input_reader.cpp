// waveinput.dat reader for CFDwavemaker -- see input_reader.h for the
// two-phase design and the robustness rules.
//
// Phase 1 collects every [section] into plain storage (inputdata plus
// passive fields of the wave objects). Phase 2 resolves the wave type,
// runs all cross-section checks and initializes the objects in a
// fixed sequence. Section order and keyword order therefore do not
// matter; unknown sections/keywords warn; malformed numbers are hard
// errors; duplicated sections are hard errors.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "input_reader.h"
#include "CFDwavemaker.h"
#include "Irregular.h"
#include "Stokes5.h"
#include "FentonStream.h"
#include "Wavemaker.h"
#include "lsgrid_spline.h"
#include "sigma_fft.h"
#include "mpi_utils.h"
#include "Utils.h"

#if defined(SWD_enable)
#include "SpectralWaveData.h"
#endif
#if defined(VTK_enable)
#include "VTKreader.h"
#endif

#include <new>       // placement new (singleton reset on re-initialization)

// application singletons (defined in CFDwavemaker.cpp)
extern CFDwavemakerInputdata inputdata;
extern Irregular irregular;
extern Stokes5 stokes5;
extern FentonStream fenton;
extern Wavemaker wavemaker;
extern lsGridSpline sgrids;
extern Ramp ramp;
#if defined(SWD_enable)
extern SpectralWaveData* swd;
#endif
#if defined(VTK_enable)
extern VTKreader vtkreader;
#endif

namespace {

// ------------------------------------------------------------ helpers

/* A sorting function for vectors which returns the indices after sorting */
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& v) {
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);
	// stable_sort keeps the file order of equal elements (several
	// directions sharing one frequency must stay paired with their
	// phases deterministically)
	std::stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
	return idx;
}

inline void ltrim(std::string& s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		return !isspace(ch);
		}));
}
inline void rtrim(std::string& s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		return !isspace(ch);
		}).base(), s.end());
}
inline void trim(std::string& s) {
	ltrim(s);
	rtrim(s);
}

// Hard error. Printed on *every* rank (not through cfdwm_mpi::err()): the input
// file is normally identical on all ranks so this is usually N copies of the
// same line, but a rank-local failure (unreadable file on one node) must not be
// swallowed. Prefixed with the rank so the copies are attributable.
[[noreturn]] void fatal(const std::string& msg) {
	if (cfdwm_mpi::size() > 1)
		std::cerr << "[rank " << cfdwm_mpi::rank() << "] ";
	std::cerr << msg << std::endl;
	cfdwm_mpi::abort(-1);
	exit(-1); // not reached
}

// first whitespace-separated token of a line ("" if none)
std::string first_token(const std::string& line) {
	std::istringstream ss(line);
	std::string tok;
	ss >> tok;
	return tok;
}

// does the line look like a numeric table row?
bool numeric_start(const std::string& line) {
	if (line.empty()) return false;
	char c = line[0];
	return isdigit((unsigned char)c) || c == '-' || c == '+' || c == '.';
}

// Strict scalar extraction: everything after the keyword must parse as
// the requested type. Malformed numbers are a hard error (they used to
// become silent zeros).
template <typename T>
T parse_value(const std::string& line, const char* section) {
	std::istringstream ss(line);
	std::string key;
	T val;
	ss >> key;
	if (!(ss >> val))
		fatal("INPUTFILE ERROR: could not parse a value for keyword '"
			+ key + "' in " + section + " (line: \"" + line + "\").");
	return val;
}

// n strict doubles from a line (keyword first)
void parse_values(const std::string& line, const char* section,
	double* out, int n) {
	std::istringstream ss(line);
	std::string key;
	ss >> key;
	for (int i = 0; i < n; i++) {
		if (!(ss >> out[i]))
			fatal(std::string("INPUTFILE ERROR: keyword '") + key + "' in "
				+ section + " needs " + std::to_string(n)
				+ " numeric values (line: \"" + line + "\").");
	}
}

// n strict doubles from a table row (no keyword). Extra trailing
// columns are allowed (and ignored), as before.
void parse_row(const std::string& line, const char* section,
	double* out, int n) {
	std::istringstream ss(line);
	for (int i = 0; i < n; i++) {
		if (!(ss >> out[i]))
			fatal(std::string("INPUTFILE ERROR: malformed table row in ")
				+ section + " - expected " + std::to_string(n)
				+ " numeric columns (line: \"" + line + "\").");
	}
}

// number of numeric columns in a table row (up to a '#' comment)
int count_columns(const std::string& line) {
	std::istringstream ss(line);
	std::string tok;
	int n = 0;
	while (ss >> tok) {
		if (tok[0] == '#') break;
		n++;
	}
	return n;
}

// string value after the keyword
std::string parse_string(const std::string& line, const char* section) {
	std::istringstream ss(line);
	std::string key, val;
	ss >> key;
	if (!(ss >> val))
		fatal("INPUTFILE ERROR: keyword '" + key + "' in " + section
			+ " needs a value.");
	return val;
}

void warn_unknown(const std::string& section, const std::string& line) {
	cfdwm_mpi::err() << "WARNING: unknown keyword '" << first_token(line)
		<< "' in " << section << " - ignored. Check the spelling "
		"against the manual." << std::endl;
}

void require_unique(bool& flag, const char* section) {
	if (flag)
		fatal(std::string("INPUTFILE ERROR: duplicate ") + section
			+ " section - each section may only appear once.");
	flag = true;
}

// one [section] with its body lines (comments/blank lines are already
// stripped by the file preprocessing)
struct Section {
	std::string name;
	std::vector<std::string> lines;
};

// -------------------------------------------------- phase 1: parsing

void parse_wave_type(const Section& s) {
	require_unique(inputdata.wavetype_read, "[wave type]");
	cfdwm_mpi::out() << "----------" << std::endl;
	cfdwm_mpi::out() << "Wave type:" << std::endl;
	cfdwm_mpi::out() << "----------" << std::endl;
	if (s.lines.empty())
		fatal("INPUTFILE ERROR: [wave type] section is empty.");
	const std::string& lineA = s.lines[0];
	if (!lineA.compare("irregular") || !lineA.compare("1")) {
		inputdata.wavetype_base = 1;
		cfdwm_mpi::out() << "Irregular perturbation wave theory specified" << std::endl;
	}
	else if (!lineA.compare("piston")) {
		inputdata.wavetype_base = 11;
		cfdwm_mpi::out() << "Wave maker theory specified" << std::endl;
	}
	else if (!lineA.compare("regular")) {
		inputdata.wavetype_base = 21;
		cfdwm_mpi::out() << "Regular 5th order Stokes wave specified" << std::endl;
	}
	else if (!lineA.compare("swd")) {
		inputdata.wavetype_base = 31;
		cfdwm_mpi::out() << "Spectral wave data (swd) specified" << std::endl;
	}
	else if (!lineA.compare("vtk")) {
		inputdata.wavetype_base = 41;
		cfdwm_mpi::out() << "Special branch VTK file interpolation" << std::endl;
	}
	else if (!lineA.compare("slope")) {
		inputdata.wavetype_base = 51;
		cfdwm_mpi::out() << "Slope of water. Water at rest" << std::endl;
	}
	else {
		fatal("INPUTFILE ERROR: Unknown wave type specified. Valid "
			"alternatives are: irregular, regular, piston, swd, vtk, slope.");
	}
	for (size_t i = 1; i < s.lines.size(); i++)
		warn_unknown(s.name, s.lines[i]);
}

void parse_general(const Section& s) {
	require_unique(inputdata.general_read, "[general input data]");
	cfdwm_mpi::out() << "-------------------" << std::endl;
	cfdwm_mpi::out() << "General input data:" << std::endl;
	cfdwm_mpi::out() << "-------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "depth") {
			inputdata.depth = parse_value<double>(lineA, s.name.c_str());
			inputdata.depth_read = true;
			cfdwm_mpi::out() << "Water depth: " << inputdata.depth << "m" << std::endl;
		}
		else if (key == "mtheta") {
			inputdata.mtheta = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Mean wave direction: " << inputdata.mtheta << " degrees" << std::endl;
		}
		else if (key == "normalize") {
			irregular.normalize = (parse_value<int>(lineA, s.name.c_str()) != 0);
			cfdwm_mpi::out() << "Normalize wave amplitudes by spectral zeroth moment: " << irregular.normalize << std::endl;
		}
		else if (key == "amplify") {
			inputdata.ampl = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Amplify (gain): " << inputdata.ampl << "m" << std::endl;
		}
		else if (key == "swl") {
			inputdata.swl = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Still water line set to: " << inputdata.swl << "m" << std::endl;
		}
		else if (key == "gravity") {
			inputdata.gravity = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Gravity set to: " << inputdata.gravity << "m/s^2" << std::endl;
		}
		else if (key == "rho") {
			inputdata.rho = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Water density rho set to: " << inputdata.rho << "kg/m^3" << std::endl;
		}
		else if (key == "nhpres_include_kinetic") {
			std::string v = parse_string(lineA, s.name.c_str());
			if (!v.compare(0, 4, "true")) {
				inputdata.nhpres_include_kinetic = true;
				cfdwm_mpi::out() << "Nonhydrostatic pressure: full Bernoulli form (incl. kinetic energy term)" << std::endl;
			}
			else {
				inputdata.nhpres_include_kinetic = false;
				cfdwm_mpi::out() << "Nonhydrostatic pressure: linearized form (kinetic energy term excluded)" << std::endl;
			}
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_refpoint(const Section& s) {
	require_unique(inputdata.refpoint_read, "[wave reference point]");
	cfdwm_mpi::out() << "---------------------" << std::endl;
	cfdwm_mpi::out() << "Wave reference point:" << std::endl;
	cfdwm_mpi::out() << "---------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "time") {
			inputdata.tofmax = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "t0: " << inputdata.tofmax << " sec" << std::endl;
		}
		else if (key == "x") {
			inputdata.x_pos = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "x0: " << inputdata.x_pos << " m" << std::endl;
		}
		else if (key == "y") {
			inputdata.y_pos = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "y0: " << inputdata.y_pos << " m" << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_ramps(const Section& s) {
	require_unique(inputdata.ramps_read, "[ramps]");
	cfdwm_mpi::out() << "------" << std::endl;
	cfdwm_mpi::out() << "Ramps:" << std::endl;
	cfdwm_mpi::out() << "------" << std::endl;
	ramp.ramp_init = true;
	double v[2];
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "time_rampup") {
			parse_values(lineA, s.name.c_str(), v, 2);
			ramp.ramp_init_time_up = true;
			ramp.time_rampup_start = v[0];
			ramp.time_rampup_end = v[1];
			cfdwm_mpi::out() << "time-rampup -- starttime: " << ramp.time_rampup_start << " sec, endtime: " << ramp.time_rampup_end << " sec." << std::endl;
		}
		else if (key == "time_rampdown") {
			parse_values(lineA, s.name.c_str(), v, 2);
			ramp.ramp_init_time_down = true;
			ramp.time_rampdown_start = v[0];
			ramp.time_rampdown_end = v[1];
			cfdwm_mpi::out() << "time-rampdown -- starttime: " << ramp.time_rampdown_start << " sec, endtime: " << ramp.time_rampdown_end << " sec." << std::endl;
		}
		else if (key == "x_rampup") {
			parse_values(lineA, s.name.c_str(), v, 2);
			ramp.ramp_init_x_up = true;
			ramp.x_rampup_start = v[0];
			ramp.x_rampup_end = v[1];
			cfdwm_mpi::out() << "x-rampup -- startpos: " << ramp.x_rampup_start << " m, endpos: " << ramp.x_rampup_end << " m." << std::endl;
		}
		else if (key == "x_rampdown") {
			parse_values(lineA, s.name.c_str(), v, 2);
			ramp.ramp_init_x_down = true;
			ramp.x_rampdown_start = v[0];
			ramp.x_rampdown_end = v[1];
			cfdwm_mpi::out() << "x-rampdown -- startpos: " << ramp.x_rampdown_start << " m, endpos: " << ramp.x_rampdown_end << " m." << std::endl;
		}
		else if (key == "y_rampup") {
			parse_values(lineA, s.name.c_str(), v, 2);
			ramp.ramp_init_y_up = true;
			ramp.y_rampup_start = v[0];
			ramp.y_rampup_end = v[1];
			cfdwm_mpi::out() << "y-rampup -- startpos: " << ramp.y_rampup_start << " m, endpos: " << ramp.y_rampup_end << " m." << std::endl;
		}
		else if (key == "y_rampdown") {
			parse_values(lineA, s.name.c_str(), v, 2);
			ramp.ramp_init_y_down = true;
			ramp.y_rampdown_start = v[0];
			ramp.y_rampdown_end = v[1];
			cfdwm_mpi::out() << "y-rampdown -- startpos: " << ramp.y_rampdown_start << " m, endpos: " << ramp.y_rampdown_end << " m." << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_irregular(const Section& s) {
	require_unique(inputdata.irregular_read, "[wave input irregular]");
	cfdwm_mpi::out() << "--------------------------" << std::endl;
	cfdwm_mpi::out() << "Irregular input data:" << std::endl;
	cfdwm_mpi::out() << "--------------------------" << std::endl;

	// scan: option keywords may appear anywhere in the section; the
	// numeric rows (frequency table, then direction table for
	// ndir > 0) are collected in file order
	bool nfreq_read = false, ndir_read = false;
	std::vector<const std::string*> rows;
	for (const std::string& lineA : s.lines) {
		if (numeric_start(lineA)) {
			rows.push_back(&lineA);
			continue;
		}
		const std::string key = first_token(lineA);
		if (key == "nfreq") {
			irregular.nfreq = parse_value<int>(lineA, s.name.c_str());
			nfreq_read = true;
		}
		else if (key == "ndir") {
			// double summation (separate direction table) is removed;
			// "ndir 0" is accepted for backward compatibility
			const int nd = parse_value<int>(lineA, s.name.c_str());
			if (nd > 0)
				fatal("INPUTFILE ERROR: ndir > 0 (double summation with a "
					"separate direction table) has been removed. Use the "
					"single-summation format: one row per component with "
					"its own direction column (WaveForge writes this "
					"format).");
			ndir_read = true;
		}
		else if (key == "kinematics_surf_extrap") {
			irregular.extrapolation_met = parse_value<int>(lineA, s.name.c_str());
			if (irregular.extrapolation_met == 2) {
				irregular.order = 2;
				cfdwm_mpi::out() << "Second order kinematics used with Taylor expansion above SWL. " << std::endl;
			}
			else if (irregular.extrapolation_met == 1) {
				irregular.order = 0;
				cfdwm_mpi::out() << "Linear wave theory used with constant value above swl." << std::endl;
			}
			else if (irregular.extrapolation_met == 0) {
				irregular.order = 0;
				cfdwm_mpi::out() << "Linear wave theory used with exponential extrapolation above z=swl." << std::endl;
			}
			else {
				cfdwm_mpi::out() << "Unknown extrapolation method. try a value between 0 and 2." << std::endl;
				cfdwm_mpi::abort(-1);
			}
		}
		else if (key == "cutoff") {
			irregular.dw_cutoff = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Second order high-frequency cut-off set to: " << irregular.dw_cutoff << std::endl;
		}
		else if (key == "equal_dk") {
			irregular.equal_dk = (parse_value<int>(lineA, s.name.c_str()) != 0);
			if (irregular.equal_dk) {
				cfdwm_mpi::out() << "Components declared on an equal-dk wavenumber grid: "
					<< "factorized second-order evaluation will be used if applicable." << std::endl;
			}
		}
		else if (key == "bandwidth") {
			std::string v = parse_string(lineA, s.name.c_str());
			if (!v.compare(0, 3, "off")) {
				cfdwm_mpi::out() << "Bandwidth: off" << std::endl;
			}
			else if (!v.compare(0, 4, "auto")) {
				cfdwm_mpi::out() << "Bandwidth: auto" << std::endl;
				inputdata.bw_auto_calc = true;
			}
			else {
				irregular.dw_bandwidth = atof(v.c_str());
				cfdwm_mpi::out() << "Bandwidth: " << irregular.dw_bandwidth << " rad/s" << std::endl;
			}
		}
		else if (key == "individual_ramp") {
			std::string v = parse_string(lineA, s.name.c_str());
			if (!v.compare(0, 4, "true")) {
				cfdwm_mpi::out() << "Individual rampup of each frequency component set to true. Three additional columns of data required when specifying frequency components." << std::endl;
				irregular.individual_ramp = true;
			}
		}
		else if (key == "vertical_theory") {
			std::string v = parse_string(lineA, s.name.c_str());
			inputdata.vertical_theory_set = true;
			if (!v.compare(0, 4, "smit")) {
				irregular.vertical_theory = 1;
				cfdwm_mpi::out() << "Smit et al. (2017) s-coordinate vertical theory selected "
					"(sigma-transformed second order; validated in v3.2.0 by the "
					"epsilon-scaling program, see sigma_developments/)." << std::endl;
			}
			else {
				irregular.vertical_theory = 0;
				cfdwm_mpi::out() << "Conventional Sharma & Dean vertical theory (default)." << std::endl;
			}
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}

	if (!nfreq_read)
		fatal("INPUTFILE ERROR: parameter nfreq missing or not specified correctly");
	(void)ndir_read;   // ndir is deprecated (single summation only)

	const int nfreq = irregular.nfreq;
	irregular.ndir = 1;

	// Single-summation component table: one row per component with its
	// own direction. The column count of the first row decides the
	// format:
	//   new    (4 cols): OMEGA A PHASE THETA -- or, with equal_dk 1,
	//                    K A PHASE THETA (the lattice is exact in k;
	//                    omega follows from the dispersion relation)
	//   legacy (5 cols): OMEGA A K PHASE THETA (k is overspecified and
	//                    recomputed from omega; deprecated)
	// individual_ramp appends 3 ramp columns to either format.
	cfdwm_mpi::out() << "Irregular seas, one directional component for each frequency specified" << std::endl;
	cfdwm_mpi::out() << "Number of frequency components: " << nfreq << std::endl;
	if ((int)rows.size() != nfreq)
		fatal("INPUTFILE ERROR: expected " + std::to_string(nfreq)
			+ " frequency component rows in [wave input irregular], "
			"found " + std::to_string(rows.size())
			+ ". Please check to make sure that the number of lines "
			"matches the specified number of frequency components");

	const int ramp_cols = irregular.individual_ramp ? 3 : 0;
	const int found_cols = count_columns(*rows[0]);
	bool legacy5;
	if (found_cols == 4 + ramp_cols)
		legacy5 = false;
	else if (found_cols == 5 + ramp_cols)
		legacy5 = true;
	else
		fatal("INPUTFILE ERROR: [wave input irregular] table rows have "
			+ std::to_string(found_cols) + " columns; expected "
			+ std::to_string(4 + ramp_cols) + " ("
			+ (irregular.equal_dk ? "K" : "OMEGA")
			+ " A PHASE THETA" + (ramp_cols ? " + 3 ramp columns" : "")
			+ ") or " + std::to_string(5 + ramp_cols)
			+ " (legacy OMEGA A K PHASE THETA).");
	const int ncols = found_cols;
	irregular.k_first = !legacy5 && irregular.equal_dk;
	irregular.k_column_read = legacy5;
	if (legacy5)
		cfdwm_mpi::out() << "NOTE: 5-column table (OMEGA A K PHASE THETA) is "
			"deprecated; the K column is recomputed from the dispersion "
			"relation. New format: " << (irregular.equal_dk ? "K" : "OMEGA")
			<< " A PHASE THETA." << std::endl;
	if (irregular.individual_ramp)
		cfdwm_mpi::out() << "Rampup start time for each frequency component "
			"specified in the trailing three columns." << std::endl;

	std::vector<double> omega, A, k, theta, phase;
	std::vector<double> f_ramptime, f_starttime, f_duration;
	double c[8];
	if (irregular.k_first)
		cfdwm_mpi::out() << "# K[1/m]    A[m]      Phase[rad]    Theta[rad]" << std::endl;
	else if (!legacy5)
		cfdwm_mpi::out() << "# OMEGA[rad/s]    A[m]      Phase[rad]    Theta[rad]" << std::endl;
	else
		cfdwm_mpi::out() << "# OMEGA[rad / s]    A[m]           K             Phase[rad]        Theta [rad]" << std::endl;
	for (int i = 0; i < nfreq; i++) {
		cfdwm_mpi::out() << *rows[i] << std::endl;
		parse_row(*rows[i], s.name.c_str(), c, ncols);
		int j = 0;
		if (legacy5) {
			omega.push_back(c[0]);
			A.push_back(c[1]);
			k.push_back(c[2]);
			phase.push_back(c[3]);
			theta.push_back(c[4]);
			j = 5;
		}
		else {
			// first column is k (equal_dk) or omega; the other follows
			// from the dispersion relation once the depth is known
			omega.push_back(irregular.k_first ? 0.0 : c[0]);
			k.push_back(irregular.k_first ? c[0] : 0.0);
			A.push_back(c[1]);
			phase.push_back(c[2]);
			theta.push_back(c[3]);
			j = 4;
		}
		if (irregular.individual_ramp) {
			f_ramptime.push_back(c[j]);
			f_starttime.push_back(c[j + 1]);
			f_duration.push_back(c[j + 2]);
		}
	}
	// Sort ascending in the primary column (omega, or k for k_first;
	// the ordering is the same, the dispersion relation is monotone)
	for (auto i : sort_indices(irregular.k_first ? k : omega)) {
		irregular.omega.push_back(omega[i]);
		irregular.A.push_back(A[i]);
		irregular.k.push_back(k[i]);
		irregular.phase.push_back(phase[i]);
		irregular.theta.push_back(theta[i]);
		if (irregular.individual_ramp) {
			irregular.freq_ramptime.push_back(f_ramptime[i]);
			irregular.freq_starttime.push_back(f_starttime[i]);
			irregular.freq_duration.push_back(f_duration[i]);
		}
		else {
			irregular.freq_ramptime.push_back(0.);
			irregular.freq_duration.push_back(HUGE);
			irregular.freq_starttime.push_back(-SMALL);
		}
	}
	irregular.initialized = true;
	inputdata.property_read = true;
}

void parse_piston(const Section& s) {
	require_unique(inputdata.piston_read, "[wave input piston]");
	if (s.lines.empty())
		fatal("INPUTFILE ERROR: [wave input piston] section is empty.");
	{
		std::istringstream ss(s.lines[0]);
		if (!(ss >> wavemaker.n_timesteps))
			fatal("INPUTFILE ERROR: first line of [wave input piston] must "
				"be the number of timesteps.");
	}
	cfdwm_mpi::out() << "Number of timesteps: " << wavemaker.n_timesteps << std::endl;
	if ((int)s.lines.size() < wavemaker.n_timesteps + 1)
		fatal("INPUTFILE ERROR: [wave input piston] declares "
			+ std::to_string(wavemaker.n_timesteps) + " timesteps but only "
			+ std::to_string(s.lines.size() - 1) + " rows follow.");

	wavemaker.PD_time = new double[wavemaker.n_timesteps];
	wavemaker.PD_ampl = new double[wavemaker.n_timesteps];
	wavemaker.PD_velo = new double[wavemaker.n_timesteps];
	wavemaker.PD_eta = new double[wavemaker.n_timesteps];
	wavemaker.PD_velo_grad = new double[wavemaker.n_timesteps];

	double c[4];
	for (int i = 0; i < wavemaker.n_timesteps; i++) {
		parse_row(s.lines[i + 1], s.name.c_str(), c, 4);
		wavemaker.PD_time[i] = c[0];
		wavemaker.PD_ampl[i] = c[1];
		wavemaker.PD_velo[i] = c[2];
		wavemaker.PD_eta[i] = c[3];
	}
	wavemaker.gradient1D(wavemaker.PD_time, wavemaker.PD_velo, wavemaker.PD_velo_grad);
	inputdata.property_read = true;
	wavemaker.initialized = true;
	cfdwm_mpi::out() << "piston timeseries read successfully." << std::endl;
}

void parse_sloping_water(const Section& s) {
	require_unique(inputdata.sloping_read, "[sloping water]");
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "slope_direction") {
			inputdata.slope_direction = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Direction of slope (degrees):  " << inputdata.slope_direction << std::endl;
		}
		else if (key == "slope_angle") {
			inputdata.slope_angle = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Slope angle (degrees) positive value gives upward slope:  " << inputdata.slope_angle << std::endl;
		}
		else if (key == "limiting_radius") {
			inputdata.slope_radius_limiter = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Limiting radius of slope (flat outside this radius):  " << inputdata.slope_radius_limiter << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_regular(const Section& s) {
	require_unique(inputdata.regular_read, "[wave input regular]");
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "theory_type") {
			std::string tt = parse_string(lineA, s.name.c_str());
			if (!tt.compare("stream") || !tt.compare("fenton")) {
				inputdata.regular_theory = 1;
				cfdwm_mpi::out() << "Regular wave theory: Fenton stream function" << std::endl;
			}
			else {
				inputdata.regular_theory = 0;
				cfdwm_mpi::out() << "Regular wave theory: Stokes 5th order" << std::endl;
			}
		}
		else if (key == "wave_order") {
			inputdata.wave_order = parse_value<int>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Stream function Fourier order: " << inputdata.wave_order << std::endl;
		}
		else if (key == "wave_length") {
			stokes5.wave_length = parse_value<double>(lineA, s.name.c_str());
			inputdata.regular_specify_length = true;
			inputdata.regular_length_given = true;
			cfdwm_mpi::out() << "Wave length:  " << stokes5.wave_length << std::endl;
		}
		else if (key == "wave_period") {
			stokes5.wave_period = parse_value<double>(lineA, s.name.c_str());
			inputdata.regular_specify_length = false;
			inputdata.regular_period_given = true;
			cfdwm_mpi::out() << "Wave period: " << stokes5.wave_period << std::endl;
		}
		else if (key == "wave_height") {
			stokes5.wave_height = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Wave height:  " << stokes5.wave_height << std::endl;
		}
		else if (key == "current_speed") {
			stokes5.current = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "current speed:  " << stokes5.current << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
	inputdata.property_read = true;
}

void parse_swd_input(const Section& s) {
	require_unique(inputdata.swd_input_read, "[wave input swd]");
	cfdwm_mpi::out() << "-------------------------------------" << std::endl;
	cfdwm_mpi::out() << "Spectral wave data:" << std::endl;
	cfdwm_mpi::out() << "-------------------------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "swdfile") {
			inputdata.swdFileName = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Swd file specified:  " << inputdata.swdFileName << std::endl;
		}
		else if (key == "rho") {
			inputdata.rho = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "rho:  " << inputdata.rho << " kg/m^3." << std::endl;
		}
		else if (key == "nsumx") {
			inputdata.nsumx = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "nsumx:  " << inputdata.nsumx << std::endl;
		}
		else if (key == "nsumy") {
			inputdata.nsumy = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "nsumy:  " << inputdata.nsumy << std::endl;
		}
		else if (key == "impl") {
			inputdata.impl = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "impl:  " << inputdata.impl << std::endl;
		}
		else if (key == "ipol") {
			inputdata.ipol = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "ipol:  " << inputdata.ipol << std::endl;
		}
		else if (key == "norder") {
			inputdata.norder = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "norder:  " << inputdata.norder << std::endl;
		}
		else if (key == "dc_bias") {
			std::string v = parse_string(lineA, s.name.c_str());
			if (!v.compare(0, 5, "false")) {
				inputdata.dc_bias = false;
			}
			else if (!v.compare(0, 4, "true")) {
				inputdata.dc_bias = true;
			}
			else {
				inputdata.dc_bias = atof(v.c_str());
			}
			cfdwm_mpi::out() << "dc_bias:   " << inputdata.dc_bias << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_seabed(const Section& s) {
	require_unique(inputdata.seabed_read, "[seabed]");
	cfdwm_mpi::out() << "--------------------------------" << std::endl;
	cfdwm_mpi::out() << "Seabed (variable bathymetry):" << std::endl;
	cfdwm_mpi::out() << "--------------------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "file") {
			inputdata.seabed_file = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "seabed surface file: "
				<< inputdata.seabed_file << std::endl;
		}
		else if (key == "zshift") {
			inputdata.seabed_zshift = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "seabed zshift: " << inputdata.seabed_zshift
				<< " m" << std::endl;
		}
		else if (key == "xshift") {
			inputdata.seabed_xshift = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "seabed xshift: " << inputdata.seabed_xshift
				<< " m" << std::endl;
		}
		else if (key == "yshift") {
			inputdata.seabed_yshift = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "seabed yshift: " << inputdata.seabed_yshift
				<< " m" << std::endl;
		}
		else if (key == "min_depth") {
			inputdata.seabed_min_depth = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "seabed minimum depth limit: "
				<< inputdata.seabed_min_depth << " m" << std::endl;
		}
		else if (key == "kinematics") {
			inputdata.seabed_kinematics =
				(parse_value<int>(lineA, s.name.c_str()) != 0);
			cfdwm_mpi::out() << "seabed used by the wave kinematics: "
				<< (inputdata.seabed_kinematics ? "yes" : "no (geometry only)")
				<< std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_lsgrid(const Section& s) {
	require_unique(inputdata.lsgrid_read, "[lsgrid]");
	cfdwm_mpi::out() << "-----------------------------------" << std::endl;
	cfdwm_mpi::out() << "Lagrangian Stretched grid (lsgrid):" << std::endl;
	cfdwm_mpi::out() << "-----------------------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "bounds") {
			parse_values(lineA, s.name.c_str(), inputdata.lsgrid_domain, 4);
			inputdata.lsgrid_domain_given = true;
		}
		else if (key == "nx") {
			inputdata.lsgrid_nx = parse_value<int>(lineA, s.name.c_str());
			inputdata.lsgrid_nx_set = true;
		}
		else if (key == "ny") {
			inputdata.lsgrid_ny = parse_value<int>(lineA, s.name.c_str());
			inputdata.lsgrid_ny_set = true;
		}
		else if (key == "points_per_wave") {
			inputdata.lsgrid_pts_per_wave =
				parse_value<double>(lineA, s.name.c_str());
			if (inputdata.lsgrid_pts_per_wave < 2.0)
				fatal("INPUTFILE ERROR: [lsgrid] 'points_per_wave' must be at "
					"least 2 (Nyquist); 6 is the default and gives about "
					"0.5% velocity error.");
		}
		else if (key == "nl") {
			// integer, or 'auto' (also the behaviour when nl is omitted):
			// the smallest layer count whose vertical profile reconstruction
			// meets 'nl_tol' is picked at initialization from the sea state
			std::string val = parse_value<std::string>(lineA, s.name.c_str());
			std::transform(val.begin(), val.end(), val.begin(),
				[](unsigned char ch) { return (char)std::tolower(ch); });
			if (val == "auto") {
				inputdata.lsgrid_nl_set = false;
			}
			else {
				try {
					size_t used = 0;
					inputdata.lsgrid_nl = std::stoi(val, &used);
					if (used != val.size()) throw std::invalid_argument("trailing");
				}
				catch (const std::exception&) {
					fatal("INPUTFILE ERROR: [lsgrid] 'nl' expects an integer or "
						"the keyword 'auto' (line: \"" + lineA + "\").");
				}
				inputdata.lsgrid_nl_set = true;
			}
		}
		else if (key == "nl_tol") {
			// target relative L2 error of the RMS velocity profile for the
			// automatic layer count (default 2e-4)
			inputdata.lsgrid_nl_tol = parse_value<double>(lineA, s.name.c_str());
			if (!(inputdata.lsgrid_nl_tol > 0.0))
				fatal("INPUTFILE ERROR: [lsgrid] 'nl_tol' must be positive.");
		}
		else if (key == "border_left" || key == "border_right"
				|| key == "border_bottom" || key == "border_top") {
			const int which = (key == "border_left")   ? 0 :
			                  (key == "border_right")  ? 1 :
			                  (key == "border_bottom") ? 2 : 3;
			inputdata.lsgrid_border[which] =
				(parse_value<int>(lineA, s.name.c_str()) != 0);
			inputdata.lsgrid_border_mode = true;
		}
		else if (key == "t0") {
			inputdata.lsgrid_t0 = parse_value<double>(lineA, s.name.c_str());
		}
		else if (key == "dt") {
			inputdata.lsgrid_dt = parse_value<double>(lineA, s.name.c_str());
			inputdata.lsgrid_dt_set = true;
		}
		else if (key == "stretch_params") {
			fatal("INPUTFILE ERROR: the tan-based 'stretch_params a b' "
				"keyword has been replaced by the single-parameter "
				"exponential stretching 'stretch_beta <beta>', where "
				"exp(beta) is the thickness ratio between the bottom "
				"and the surface layer (0 = uniform, default 3.0). "
				"Update your [lsgrid] section.");
		}
		else if (key == "stretch_beta") {
			std::string val = parse_value<std::string>(lineA, s.name.c_str());
			std::transform(val.begin(), val.end(), val.begin(),
				[](unsigned char ch) { return (char)std::tolower(ch); });
			if (val == "auto") {
				inputdata.lsgrid_stretch_auto = true;
			}
			else {
				try {
					size_t used = 0;
					inputdata.lsgrid_stretch_beta = std::stod(val, &used);
					if (used != val.size()) throw std::invalid_argument("trailing");
				}
				catch (const std::exception&) {
					fatal("INPUTFILE ERROR: 'stretch_beta' expects a number or "
						"the keyword 'auto' (line: \"" + lineA + "\").");
				}
				inputdata.lsgrid_stretch_auto = false;
			}
		}
		else if (key == "ignore_subdomain") {
			// Removed in v3.7.6. Not silently ignored: dropping it would
			// quietly restore full-domain storage, which is a large memory
			// regression with no visible symptom.
			fatal("INPUTFILE ERROR: 'ignore_subdomain' was removed in v3.7.6. "
				"Say which borders need kinematics instead, e.g.\n"
				"    border_left 1\n    border_right 0\n"
				"    border_bottom 1\n    border_top 0\n"
				"Only a strip along the named borders is stored; the interior "
				"is dropped. Unlike ignore_subdomain this is independent of "
				"the CFD grid resolution.");
		}
		else if (key == "ignore_at_init") {
			inputdata.lsgrid_ignore_at_init =
				(parse_value<int>(lineA, s.name.c_str()) != 0);
		}
		else if (key == "store_pressure") {
			// 0/1, or 'auto' (the default): on for border-mode grids
			std::string val = parse_value<std::string>(lineA, s.name.c_str());
			std::transform(val.begin(), val.end(), val.begin(),
				[](unsigned char ch) { return (char)std::tolower(ch); });
			if (val == "auto")
				inputdata.lsgrid_store_pressure = -1;
			else if (val == "0" || val == "1")
				inputdata.lsgrid_store_pressure = (val == "1");
			else
				fatal("INPUTFILE ERROR: [lsgrid] 'store_pressure' expects 0, 1 "
					"or 'auto' (line: \"" + lineA + "\").");
		}
		else if (key == "init_only") {
			inputdata.lsgrid_init_only =
				(parse_value<int>(lineA, s.name.c_str()) != 0);
		}
		else if (key == "scheme") {
			// legacy keyword from the old linear/spline scheme choice;
			// spline interpolation is the only scheme since v3.1
			cfdwm_mpi::out() << "note: the 'scheme' keyword is deprecated - "
				"spline interpolation is always used." << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
	cfdwm_mpi::out() << "Interpolation Scheme: spline" << std::endl;
	cfdwm_mpi::out() << "LS grid domain bounds: " << std::endl;
	cfdwm_mpi::out() << "\tXMIN: " << inputdata.lsgrid_domain[0] << std::endl;
	cfdwm_mpi::out() << "\tXMAX: " << inputdata.lsgrid_domain[1] << std::endl;
	cfdwm_mpi::out() << "\tYMIN: " << inputdata.lsgrid_domain[2] << std::endl;
	cfdwm_mpi::out() << "\tYMAX: " << inputdata.lsgrid_domain[3] << std::endl;
	cfdwm_mpi::out() << "Grid resolution: " << std::endl;
	cfdwm_mpi::out() << "nx: " << inputdata.lsgrid_nx << std::endl;
	cfdwm_mpi::out() << "ny: " << inputdata.lsgrid_ny << std::endl;
	if (inputdata.lsgrid_nl_set)
		cfdwm_mpi::out() << "nl: " << inputdata.lsgrid_nl << std::endl;
	else
		cfdwm_mpi::out() << "nl: auto (smallest layer count with relative "
			"profile error <= " << inputdata.lsgrid_nl_tol
			<< "; resolved at initialization)" << std::endl;
	cfdwm_mpi::out() << "time init: " << inputdata.lsgrid_t0 << std::endl;
	cfdwm_mpi::out() << "dt: " << inputdata.lsgrid_dt << std::endl;
	if (inputdata.lsgrid_stretch_auto) {
		cfdwm_mpi::out() << "Layer stretching: auto (beta chosen at initialization to "
			"minimise the vertical interpolation error of the loaded sea state)"
			<< std::endl;
	}
	else {
		const int nlr = std::max(inputdata.lsgrid_nl, 3);
		cfdwm_mpi::out() << "Layer stretching: beta=" << inputdata.lsgrid_stretch_beta
			<< " (bottom/surface layer thickness ratio "
			<< std::exp(inputdata.lsgrid_stretch_beta * double(nlr - 2) / double(nlr - 1))
			<< " for nl=" << inputdata.lsgrid_nl << ")" << std::endl;
	}
	if (irregular.equal_dk && !inputdata.vertical_theory_set
			&& irregular.vertical_theory != 1) {
		// An equal-dk lattice is the sigma-FFT case, and the FFT fills the
		// grid with the smit sigma-transform -- Sharma & Dean would not
		// match the field that ends up stored. Requiring the keyword here
		// only creates a way to silently lose the fast path, so default it.
		irregular.vertical_theory = 1;
		cfdwm_mpi::out() << "equal_dk: vertical_theory defaults to smit (the "
			"sigma-FFT fills the grid with the smit sigma-transform; state "
			"'vertical_theory sharma' explicitly to override, which disables "
			"the FFT fast path)." << std::endl;
	}
	if (inputdata.lsgrid_border_mode) {
		cfdwm_mpi::out() << "Kinematics stored along borders:"
			<< (inputdata.lsgrid_border[0] ? " left" : "")
			<< (inputdata.lsgrid_border[1] ? " right" : "")
			<< (inputdata.lsgrid_border[2] ? " bottom" : "")
			<< (inputdata.lsgrid_border[3] ? " top" : "")
			<< " (the interior is not stored)" << std::endl;
	}
	cfdwm_mpi::out() << "Ignore subdomain at t_init: " << inputdata.lsgrid_ignore_at_init << std::endl;
}

void parse_lsgrid2vtk(const Section& s) {
	require_unique(inputdata.lsgrid2vtk_read, "[lsgrid2vtk]");
	cfdwm_mpi::out() << "--------------------" << std::endl;
	cfdwm_mpi::out() << "VTK output settings:" << std::endl;
	cfdwm_mpi::out() << "--------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "storage_path") {
			sgrids.vtk_directory_path = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Directory for storage of vtk files: " << sgrids.vtk_directory_path << std::endl;
		}
		else if (key == "filename") {
			sgrids.vtk_prefix = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "filename prefix: " << sgrids.vtk_prefix << std::endl;
		}
		else if (key == "vtk_timelabel") {
			sgrids.vtk_timelabel = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "time label to use in vtk specified as: " << sgrids.vtk_timelabel << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
}

void parse_vtk_input(const Section& s) {
	require_unique(inputdata.vtk_input_read, "[wave input vtk]");
#if defined(VTK_enable)
	cfdwm_mpi::out() << "----------------------------" << std::endl;
	cfdwm_mpi::out() << "VTK interpolation from file:" << std::endl;
	cfdwm_mpi::out() << "----------------------------" << std::endl;
	for (const std::string& lineA : s.lines) {
		const std::string key = first_token(lineA);
		if (key == "storage_path") {
			vtkreader.vtkfilepath = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Directory where vtk files are stored: " << vtkreader.vtkfilepath << std::endl;
		}
		else if (key == "filename") {
			vtkreader.vtk_prefix = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "filename prefix: " << vtkreader.vtk_prefix << std::endl;
		}
		else if (key == "name_velocity_field") {
			vtkreader.Uname = parse_string(lineA, s.name.c_str());
			cfdwm_mpi::out() << "Name of velocity scalar field: " << vtkreader.Uname << std::endl;
		}
		else if (key == "t_start") {
			vtkreader.t_start = parse_value<double>(lineA, s.name.c_str());
			cfdwm_mpi::out() << "User specified start time (t0): " << vtkreader.t_start << std::endl;
		}
		else if (key == "update_height_function_every_timestep") {
			std::string v = parse_string(lineA, s.name.c_str());
			if (!v.compare(0, 5, "false")) {
				vtkreader.recompute_betah_every_timstep = false;
			}
			else if (!v.compare(0, 4, "true")) {
				vtkreader.recompute_betah_every_timstep = true;
			}
			else {
				vtkreader.recompute_betah_every_timstep = atof(v.c_str());
			}
			cfdwm_mpi::out() << "update_height_function_every_timestep:   " << vtkreader.recompute_betah_every_timstep << std::endl;
		}
		// Multilayer vertical velocity reconstruction (only used when cell field 'h' is present):
		//   "ppr" (default, conservative parabolic) or "linear" (piecewise linear between cell centres).
		else if (key == "vertical_interp") {
			std::string v = parse_string(lineA, s.name.c_str());
			vtkreader.vertical_interp = (!v.compare(0, 6, "linear")) ? 1 : 0;
			cfdwm_mpi::out() << "Multilayer vertical velocity reconstruction: "
				<< (vtkreader.vertical_interp == 0 ? "PPR" : "linear") << std::endl;
		}
		else {
			warn_unknown(s.name, lineA);
		}
	}
#else
	(void)s;
	fatal("VTK interpolation not supported in current compiled version of CFDwavemaker. Recompile width flag -DVTK_enable=1.");
#endif
}

// split the (comment-stripped) input text into sections and dispatch.
// The @v version line has already been consumed by the caller.
void parse_all_sections(std::istringstream& f) {
	std::vector<Section> sections;
	std::string lineA;
	Section* cur = nullptr;
	while (std::getline(f, lineA)) {
		trim(lineA);
		if (lineA.empty()) continue;
		if (lineA[0] == '[') {
			sections.push_back(Section());
			cur = &sections.back();
			cur->name = lineA;
		}
		else if (cur != nullptr) {
			cur->lines.push_back(lineA);
		}
		else {
			cfdwm_mpi::err() << "WARNING: line before the first [section] tag "
				"ignored: \"" << lineA << "\"" << std::endl;
		}
	}

	for (const Section& s : sections) {
		if (s.name == "[wave type]") parse_wave_type(s);
		else if (s.name == "[general input data]") parse_general(s);
		else if (s.name == "[wave reference point]") parse_refpoint(s);
		else if (s.name == "[ramps]") parse_ramps(s);
		else if (s.name == "[wave input irregular]") parse_irregular(s);
		else if (s.name == "[wave input piston]") parse_piston(s);
		else if (s.name == "[sloping water]") parse_sloping_water(s);
		else if (s.name == "[wave input regular]") parse_regular(s);
		else if (s.name == "[wave input swd]") parse_swd_input(s);
		else if (s.name == "[seabed]") parse_seabed(s);
		else if (s.name == "[lsgrid]") parse_lsgrid(s);
		else if (s.name == "[lsgrid2vtk]") parse_lsgrid2vtk(s);
		else if (s.name == "[wave input vtk]") parse_vtk_input(s);
		else {
			cfdwm_mpi::err() << "WARNING: unknown section " << s.name
				<< " ignored (" << s.lines.size() << " lines). Known "
				"sections: [wave type], [general input data], "
				"[wave reference point], [ramps], [wave input irregular], "
				"[wave input regular], [wave input piston], "
				"[wave input swd], [wave input vtk], [sloping water], "
				"[seabed], [lsgrid], [lsgrid2vtk]. Check the spelling."
				<< std::endl;
		}
	}
}

// -------------------------------------- phase 2: resolve + initialize

/* ---------------------------------------------------------------------------
   Automatic lsgrid sizing.

   The grid only has to resolve the field it stores, and that field is fully
   described by the wave components -- so nx/ny/dt need not be stated by hand:
   the shortest wave present sets the node spacing and the shortest period sets
   the time step. LSGRID_PTS_PER_WAVE points per shortest wavelength (and per
   shortest period) is the target. The grid is read through a cubic Hermite
   interpolant, which converges at fourth order: measured against a converged
   reference, the velocity error is ~1.7% at 4 points, 0.5% at 6, 0.14% at 8.
   Six is the default -- cost grows as nx*ny, but with border storage the grid
   is a few per cent of the domain anyway, so the accuracy is worth far more
   than the memory. Override with 'points_per_wave' in [lsgrid].

   Under the sigma-FFT none of this is used -- there the lattice is fixed by
   the component set exactly (see the snapping in the [lsgrid] block). This is
   for the paths that have no lattice: direct summation, variable bathymetry,
   and SWD.
   --------------------------------------------------------------------------- */
static const double LSGRID_PTS_PER_WAVE_DEFAULT = 6.0;   // see lsgrid_pts_per_wave

/* Shortest wave (m) and shortest period (s) among the loaded components. */
static bool lsgrid_shortest_wave(const Irregular& irr, double& Lmin, double& Tmin)
{
	const int N = irr.nfreq * irr.ndir;
	if (N <= 0 || (int)irr.k.size() < N || (int)irr.omega.size() < N)
		return false;
	double kmax = 0., ommax = 0.;
	for (int i = 0; i < N; ++i) {
		kmax = std::max(kmax, irr.k[i]);
		ommax = std::max(ommax, irr.omega[i]);
	}
	if (!(kmax > 0.) || !(ommax > 0.))
		return false;
	Lmin = 2.0 * M_PI / kmax;
	Tmin = 2.0 * M_PI / ommax;
	return true;
}

/* Fill in nx/ny/dt wherever they were not given in the input file. */
static void lsgrid_autosize(double Lmin, double Tmin, const char* src)
{
	const double pts = inputdata.lsgrid_pts_per_wave;
	const double dtarget = Lmin / pts;
	const double span_x = inputdata.lsgrid_domain[1] - inputdata.lsgrid_domain[0];
	const double span_y = inputdata.lsgrid_domain[3] - inputdata.lsgrid_domain[2];
	if (!inputdata.lsgrid_nx_set)
		inputdata.lsgrid_nx = std::max(2, (int)std::ceil(span_x / dtarget) + 1);
	if (!inputdata.lsgrid_ny_set)
		// a long-crested grid keeps ny = 1; only size y when it has extent
		inputdata.lsgrid_ny = (span_y > 0.)
			? std::max(2, (int)std::ceil(span_y / dtarget) + 1) : 1;
	if (!inputdata.lsgrid_dt_set)
		inputdata.lsgrid_dt = Tmin / pts;
	if (inputdata.lsgrid_nx_set && inputdata.lsgrid_ny_set
			&& inputdata.lsgrid_dt_set)
		return;                 // nothing was derived: nothing to report
	cfdwm_mpi::out() << "lsgrid: sized from the " << src
		<< " -- shortest wave " << Lmin << " m, shortest period " << Tmin
		<< " s, at " << pts << " points each -> nx="
		<< inputdata.lsgrid_nx << ", ny=" << inputdata.lsgrid_ny << ", dt="
		<< inputdata.lsgrid_dt << " s"
		<< (inputdata.lsgrid_nx_set ? " (nx from file)" : "")
		<< (inputdata.lsgrid_ny_set ? " (ny from file)" : "")
		<< (inputdata.lsgrid_dt_set ? " (dt from file)" : "")
		<< "." << std::endl;
}

int resolve_and_initialize() {
	// ---- resolve the wave type -------------------------------------
	if (!inputdata.wavetype_read)
		fatal("INPUTFILE ERROR: no [wave type] section found.");
	inputdata.wavetype = inputdata.wavetype_base;
	if (inputdata.lsgrid_read) {
		if (inputdata.wavetype_base == 1) {
			inputdata.wavetype = 5;
			cfdwm_mpi::out() << "LS grid + Spline + irregular" << std::endl;
		}
		else if (inputdata.wavetype_base == 31) {
			inputdata.wavetype = 34;
			cfdwm_mpi::out() << "LS grid + Spline + swd" << std::endl;
		}
		else {
			fatal("INPUTFILE ERROR: LS grid may currently only be used with irregular second order wave theory or swd");
		}
	}

	// ---- cross-section consistency ---------------------------------
	if (inputdata.irregular_read && inputdata.wavetype_base != 1)
		fatal("INPUTFILE ERROR: irregular wave components does not match the specified wave type. Check inputfile");
	if (inputdata.piston_read && inputdata.wavetype_base != 11)
		fatal("piston wave maker properties does not match the specified wave type. Check inputfile");
	if (inputdata.regular_read && inputdata.wavetype_base != 21)
		fatal("INPUTFILE ERROR: Stokes wave properties does not match the specified wave type. Check inputfile");
	if (inputdata.swd_input_read && inputdata.wavetype_base != 31)
		fatal("INPUTFILE ERROR: please set [wave type] = swd when using keyword [wave input swd].");
	if (inputdata.lsgrid2vtk_read && inputdata.wavetype != 5)
		fatal("InputError: VTK output only works in combination with LSgrid at the moment.");
	if (inputdata.lsgrid2vtk_read)
		sgrids.dump_vtk = true;

	if (inputdata.wavetype == 5 && irregular.order != 2)
		fatal("LSgrid interpolation uses strictly second order theory. Doesnt make sense to use grid interpolation for linear theory since this is very fast anyway. Set kinematics_surf_extrap to 2 under [wave input irregular], or remove [lsgrid].");

	if (inputdata.seabed_read && inputdata.wavetype != 1
			&& inputdata.wavetype != 5)
		fatal("INPUTFILE ERROR: [seabed] (variable bathymetry) is "
			"only supported for irregular waves (direct evaluation or "
			"lsgrid interpolation). Other wave types assume constant "
			"depth.");
	if (inputdata.seabed_read && inputdata.depth_read
			&& inputdata.seabed_kinematics)
		fatal("INPUTFILE ERROR: 'depth' under [general input data] "
			"must not be given together with a [seabed] section - the "
			"water depth is defined by the seabed surface. Remove the "
			"depth keyword (or set 'kinematics 0' in [seabed] for "
			"geometry-only mode, where the kinematics use the given "
			"constant depth).");
	if (inputdata.seabed_read && !inputdata.seabed_kinematics
			&& !inputdata.depth_read)
		fatal("INPUTFILE ERROR: [seabed] with 'kinematics 0' "
			"(geometry-only mode) needs an explicit 'depth' under "
			"[general input data] - the wave kinematics assume that "
			"constant depth.");

	cfdwm_mpi::out() << "\n-----------------------------------------------" << std::endl;
	cfdwm_mpi::out() << "Input file read successfully." << std::endl;
	cfdwm_mpi::out() << "***********************************************\n\n" << std::endl;

	cfdwm_mpi::out() << "Wavetype: " << inputdata.wavetype << std::endl;

	// ---- irregular waves (direct or lsgrid) ------------------------
	if (inputdata.wavetype == 1 || inputdata.wavetype == 5) {
		if (inputdata.seabed_read) {
			if (inputdata.seabed_file.empty())
				fatal("INPUTFILE ERROR: [seabed] section given "
					"without a 'file <path>' keyword.");
			irregular.seabed.set_min_depth(inputdata.seabed_min_depth);
			if (!irregular.seabed.load(inputdata.seabed_file,
					inputdata.seabed_zshift, inputdata.seabed_xshift,
					inputdata.seabed_yshift))
				fatal("INPUTFILE ERROR: failed to load seabed "
					"surface " + inputdata.seabed_file);
			sgrids.seabed = &irregular.seabed;
			irregular.seabed.set_kinematics(inputdata.seabed_kinematics);
			cfdwm_mpi::out() << "Seabed surface loaded: "
				<< irregular.seabed.ntriangles() << " triangles, depth range ["
				<< irregular.seabed.min_depth() << ", "
				<< irregular.seabed.max_depth() << "] m (before the "
				<< inputdata.seabed_min_depth << " m minimum-depth clamp)."
				<< std::endl;
			if (inputdata.seabed_kinematics) {
				// representative depth (at the wave reference point) for
				// the spectral summaries and as a fallback; the
				// kinematics use the local depth everywhere
				inputdata.depth = irregular.seabed.depth(inputdata.x_pos,
					inputdata.y_pos);
				cfdwm_mpi::out() << "Reference depth at the wave reference point: "
					<< inputdata.depth << " m." << std::endl;
				cfdwm_mpi::out() << "Slope-corrected first-order kinematics active "
					"(local depth + dispersion, exact bottom BC)." << std::endl;
				// Both of these describe the ACTIVE-bathymetry kinematics, so
				// they belong inside this branch: in geometry-only mode the
				// second order is the ordinary constant-depth Sharma & Dean
				// and smit is honoured, not overridden. Printing them
				// unconditionally (as before) contradicted the geometry-only
				// line right above and the sigma-FFT log right below.
				if (irregular.order >= 2)
					cfdwm_mpi::out() << "Second order: depth-local Sharma & Dean "
						"(interaction-limited pair set, evaluated at the local "
						"depth)." << std::endl;
				if (irregular.vertical_theory == 1)
					cfdwm_mpi::out() << "WARNING: vertical_theory smit is ignored "
						"with [seabed] kinematics ON (the slope-corrected "
						"eigenfunction is used instead)." << std::endl;
			}
			else
				cfdwm_mpi::out() << "Seabed kinematics OFF (geometry-only): the "
					"wave kinematics assume the constant depth "
					<< inputdata.depth << " m from [general input data] and the "
					"requested vertical theory is honoured in full; the surface "
					"is served to linked solvers via wave_Seabed*()."
					<< std::endl;
		}
		irregular.depth = inputdata.depth;
		irregular.mtheta = inputdata.mtheta;
		irregular.tofmax = inputdata.tofmax;
		irregular.fpoint[0] = inputdata.x_pos;
		irregular.fpoint[1] = inputdata.y_pos;
		irregular.swl = inputdata.swl;
		irregular.g = inputdata.gravity;
		irregular.rho = inputdata.rho;
		irregular.ampl = inputdata.ampl;

		// One of omega/k is primary, the other follows from the
		// dispersion relation at the reference depth (known only now):
		//  - k_first (new equal_dk format): omega = sqrt(g k tanh(k h)),
		//    closed form; the k lattice stays exact as read.
		//  - otherwise omega is primary: recompute k (Newton, machine
		//    precision) and, for legacy 5-column files, warn if the file
		//    value disagreed (usually a different generation depth).
		//    With [seabed] the reference k is re-solved locally anyway.
		if (irregular.k_first) {
			for (int ci = 0; ci < (int)irregular.k.size(); ci++) {
				const double kh = irregular.k[ci] * irregular.depth;
				const double th = (kh > 20.0) ? 1.0 : std::tanh(kh);
				irregular.omega[ci] = std::sqrt(
					irregular.g * irregular.k[ci] * th);
			}
		}
		else {
			int nbad = 0;
			double worst = 0.;
			int iworst = -1;
			for (int ci = 0; ci < (int)irregular.k.size(); ci++) {
				const double k_new = irregular.dispersion_k(irregular.omega[ci],
					irregular.depth);
				const double rel = std::fabs(k_new - irregular.k[ci])
					/ std::max(k_new, 1e-30);
				if (rel > 1e-3) {
					nbad++;
					if (rel > worst) { worst = rel; iworst = ci; }
				}
				irregular.k[ci] = k_new;
			}
			if (nbad > 0 && irregular.k_column_read
					&& !irregular.seabed.active()) {
				cfdwm_mpi::out() << "WARNING: the wavenumber column disagrees "
					"with the dispersion relation at depth " << irregular.depth
					<< " m for " << nbad << " of " << irregular.k.size()
					<< " components (worst: component " << (iworst + 1)
					<< ", " << (100. * worst) << " %). The recomputed "
					"wavenumbers are used. Check that the wave input was "
					"generated for this depth." << std::endl;
			}
		}

		if (inputdata.bw_auto_calc) {
			irregular.dw_bandwidth = irregular.bandwidth_estimator();
			cfdwm_mpi::out() << "BW_autocalc: Bandwidth parameter set to " << irregular.dw_bandwidth << " rad/s." << std::endl;
		}
		irregular.normalize_data();
		irregular.calculate_bwindices();
		if (cfdwm_mpi::rank() == 0)   // MPI: one writer for the dump file
			irregular.dumpSpectralComponents();
		irregular.print_spectral_parameters();
		cfdwm_mpi::out() << "irregorder: " << irregular.order << std::endl;
		if (irregular.order == 2) {
			cfdwm_mpi::out() << "done initializing" << std::endl;
			irregular.rebuild_second_order_cache();
			cfdwm_mpi::out() << "done initializing" << std::endl;
		}

		if (!inputdata.property_read)
			fatal("INPUTFILE ERROR: Irregular wave selected, but no wave components/wave properties found in input file.");
	}
	else if (inputdata.wavetype == 11) {
		if (!inputdata.property_read)
			fatal("INPUTFILE ERROR: Wave maker selected, but no wavemaker property data found in input file.");
	}
	else if (inputdata.wavetype == 21) {
		if (!inputdata.property_read)
			fatal("INPUTFILE ERROR: Regular wave selected, but no wave properties found in input file.");
		// Length <-> period conversion, deferred from the parse phase so
		// it always sees the final depth/gravity regardless of whether
		// [general input data] came before or after [wave input regular]
		if (inputdata.regular_length_given && inputdata.regular_period_given)
			fatal("INPUTFILE ERROR: both wave_length and wave_period are "
				"given under [wave input regular]. Specify either one - "
				"the other follows from the dispersion relation.");
		stokes5.depth = inputdata.depth;
		stokes5.gravity = inputdata.gravity;
		if (inputdata.regular_length_given) {
			stokes5.wave_period = stokes5.wavePeriodFromLength(stokes5.wave_length);
			cfdwm_mpi::out() << "Wave period computed from specified wave length: " << stokes5.wave_period << std::endl;
		}
		else if (inputdata.regular_period_given) {
			stokes5.wave_length = stokes5.wavelengthFromPeriod(stokes5.wave_period);
			cfdwm_mpi::out() << "Wave length computed from specified wave period: " << stokes5.wave_length << std::endl;
		}
		stokes5.initialized = true;

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
			fenton.current = stokes5.current;   // shares the current_speed keyword

			// Solve: specify_length uses wave_length, else wave_period. Current is Eulerian mean.
			fenton.set_stream_properties(stokes5.wave_height, stokes5.wave_length, stokes5.wave_period,
				inputdata.regular_specify_length, /*current_type=Eulerian*/0,
				inputdata.wave_order, /*height-continuation steps*/5);
			cfdwm_mpi::out() << "Fenton stream function wave initialized (N=" << inputdata.wave_order
				<< ", L=" << fenton.wave_length << ", T=" << fenton.wave_period
				<< ", c=" << fenton.c << ")." << std::endl;
		}
		else {
			// Stokes 5th order
			stokes5.theta = inputdata.mtheta;
			stokes5.x0 = inputdata.x_pos;
			stokes5.y0 = inputdata.y_pos;
			stokes5.t0 = inputdata.tofmax;
			stokes5.z0 = inputdata.swl;
			stokes5.rho = inputdata.rho;

			stokes5.set_stokes5_properties(stokes5.wave_length, stokes5.wave_height);
			cfdwm_mpi::out() << "Stokes 5th wave initialized and ready to go." << std::endl;
		}
	}

	// ---- lsgrid over irregular waves -------------------------------
	if (inputdata.wavetype == 5) {
		// copy the parsed parameters onto the grid object
		if (!inputdata.lsgrid_domain_given && !inputdata.lsgrid_domain_override)
			fatal("INPUTFILE ERROR: the [lsgrid] extent is not defined. Give "
				"it either as\n    bounds <xmin> <xmax> <ymin> <ymax>\n"
				"in the [lsgrid] section, or by calling\n"
				"    wave_SetLsgridBounds(xmin, xmax, ymin, ymax)\n"
				"from the linked solver BEFORE wave_Initialize(). The latter "
				"is preferable when the solver already knows its own domain, "
				"since the two cannot then drift apart.");
		if (inputdata.lsgrid_domain_override)
			for (int b = 0; b < 4; b++)
				inputdata.lsgrid_domain[b] = inputdata.lsgrid_domain_set[b];
		{
			double Lmin, Tmin;
			if (lsgrid_shortest_wave(irregular, Lmin, Tmin))
				lsgrid_autosize(Lmin, Tmin, "wave component set");
		}
		if (inputdata.lsgrid_border_override) {
			for (int b = 0; b < 4; b++)
				inputdata.lsgrid_border[b] = inputdata.lsgrid_border_set[b];
			inputdata.lsgrid_border_mode = true;
		}
		sgrids.domain[0] = inputdata.lsgrid_domain[0];
		sgrids.domain[1] = inputdata.lsgrid_domain[1];
		sgrids.domain[2] = inputdata.lsgrid_domain[2];
		sgrids.domain[3] = inputdata.lsgrid_domain[3];
		sgrids.nx = inputdata.lsgrid_nx;
		sgrids.ny = inputdata.lsgrid_ny;
		sgrids.nl = inputdata.lsgrid_nl;
		sgrids.t0 = inputdata.lsgrid_t0;
		sgrids.dt = inputdata.lsgrid_dt;
		sgrids.stretch_beta = inputdata.lsgrid_stretch_beta;
		sgrids.stretch_auto = inputdata.lsgrid_stretch_auto;
		sgrids.nl_auto = !inputdata.lsgrid_nl_set;
		sgrids.nl_tol = inputdata.lsgrid_nl_tol;
		for (int b = 0; b < 4; b++) sgrids.border[b] = inputdata.lsgrid_border[b];
		sgrids.border_margin = inputdata.lsgrid_border_margin;
		sgrids.border_mode = inputdata.lsgrid_border_mode;
		sgrids.init_only = inputdata.lsgrid_init_only;
		sgrids.ignore_at_init = inputdata.lsgrid_ignore_at_init;
		// Pressure field: 'store_pressure' explicit, else auto -- ON for
		// border-mode grids (a GABC subtracts the incident non-hydrostatic
		// pressure there, and without a stored field wave_nhPres() can only
		// serve a first-order fallback), OFF for full-domain initialization
		// grids (initialization consumes eta/u/v/w only, and the field
		// costs ~1/3 more memory and fill time).
		sgrids.store_pressure = (inputdata.lsgrid_store_pressure >= 0)
			? (inputdata.lsgrid_store_pressure != 0)
			: inputdata.lsgrid_border_mode;
		if (sgrids.store_pressure)
			cfdwm_mpi::out() << "lsgrid: the dynamic pressure is stored as a "
				"grid field (store_pressure "
				<< (inputdata.lsgrid_store_pressure >= 0
					? "set" : "auto: border mode")
				<< "); wave_nhPres()/wave_DynPres() are served by spline "
				"interpolation at full second order." << std::endl;
		sgrids.allocate();

		sgrids.water_depth = inputdata.depth;
		sgrids.swl = inputdata.swl;

		// equal-dk + smit: lock the lsgrid to the spatial period(s) of the
		// wave-vector lattice. The field repeats with Lx = 2*pi/dkx (and
		// Ly = 2*pi/dky for 2D short-crested lattices), so the grid
		// resolves the components exactly and queries outside the domain
		// wrap around periodically.
		{
			bool twod = false;
			double dkx = 0.0, dky = 0.0;
			int jxmax = 0, jymax = 0;

			// A 2D wave-vector lattice defines the directions absolutely:
			// mtheta would rotate it off the grid axes and break the FFT
			// structure. If removing mtheta reveals a 2D lattice, warn and
			// ignore the value (the generator folds the mean direction
			// into the component angles).
			if (irregular.equal_dk && irregular.vertical_theory == 1
				&& irregular.mtheta != 0.0 && !irregular.seabed.active()) {
				const double mtheta_saved = irregular.mtheta;
				irregular.mtheta = 0.0;
				irregular.rebuild_second_order_cache();
				if (SigmaFFTFill::probe(irregular, inputdata.swl, twod, dkx,
						dky, jxmax, jymax) && twod) {
					cfdwm_mpi::err() << "WARNING: mtheta = " << mtheta_saved
						<< " deg is IGNORED: the components form a 2D "
						"wave-vector lattice, which defines the wave "
						"directions absolutely. Fold the mean direction "
						"into the component directions in the generator "
						"(WaveForge does this automatically)."
						<< std::endl;
					inputdata.mtheta = 0.0;
				}
				else {
					irregular.mtheta = mtheta_saved;
					irregular.rebuild_second_order_cache();
				}
			}

			if (!irregular.seabed.active()
					&& SigmaFFTFill::probe(irregular, inputdata.swl, twod, dkx,
						dky, jxmax, jymax)) {
				// Snap the lsgrid to the FFT nodes of the equal-dk
				// lattice. The FFT evaluates one full spatial period
				// anchored at the grid origin, but only the nodes inside
				// the requested window are stored: the grid spacing is
				// pinned to the tile node spacing (every stored value is
				// an exact collocation value), xmin/ymin stay as given,
				// and xmax/ymax round outward to the next node. The
				// file's nx/ny act as a lower bound on the per-period
				// FFT resolution. A window covering the full period
				// locks to exactly one tile and is periodic.
				const double Lpx = 2.0 * PI / dkx;
				int nfftx = 1;
				const int floorx = inputdata.lsgrid_nx_set
					? std::max(2 * (2 * jxmax) + 2, sgrids.nx - 1)
					: 2 * (2 * jxmax) + 2;   // components alone
				while (nfftx < floorx) nfftx <<= 1;
				const double dxs = Lpx / nfftx;
				const double span_x = sgrids.domain[1] - sgrids.domain[0];
				if (span_x >= Lpx - 0.5 * dxs) {
					sgrids.nx = nfftx + 1;
					sgrids.domain[1] = sgrids.domain[0] + Lpx;
					sgrids.periodic_x = true;
				}
				else {
					const int nwx = std::max(2, (int)std::ceil(
						span_x / dxs - 1.0e-9) + 1);
					sgrids.nx = nwx;
					sgrids.domain[1] = sgrids.domain[0] + (nwx - 1) * dxs;
					sgrids.periodic_x = false;
				}
				if (twod) {
					const double Lpy = 2.0 * PI / dky;
					int nffty = 1;
					const int floory = inputdata.lsgrid_ny_set
						? std::max(2 * (2 * jymax) + 2, sgrids.ny - 1)
						: 2 * (2 * jymax) + 2;
					while (nffty < floory) nffty <<= 1;
					const double dys = Lpy / nffty;
					const double span_y =
						sgrids.domain[3] - sgrids.domain[2];
					if (span_y >= Lpy - 0.5 * dys) {
						sgrids.ny = nffty + 1;
						sgrids.domain[3] = sgrids.domain[2] + Lpy;
						sgrids.periodic_y = true;
					}
					else {
						const int nwy = std::max(2, (int)std::ceil(
							span_y / dys - 1.0e-9) + 1);
						sgrids.ny = nwy;
						sgrids.domain[3] =
							sgrids.domain[2] + (nwy - 1) * dys;
						sgrids.periodic_y = false;
					}
				}
				sgrids.reallocate();
				cfdwm_mpi::out() << "sigma-fft: lsgrid snapped to the equal-dk "
					"lattice nodes: x in [" << sgrids.domain[0] << ", "
					<< sgrids.domain[1] << "] nx=" << sgrids.nx
					<< (sgrids.periodic_x ? " (full period, periodic)"
						: " (window)");
				if (twod) {
					cfdwm_mpi::out() << ", y in [" << sgrids.domain[2] << ", "
						<< sgrids.domain[3] << "] ny=" << sgrids.ny
						<< (sgrids.periodic_y ? " (full period, periodic)"
							: " (window)");
				}
				cfdwm_mpi::out() << ", node spacing dx=" << dxs << " m."
					<< std::endl;
			}
		}

		sgrids.set_ignore();

		if (sgrids.ignore_at_init) {
			// The ignored region is never filled, so it need never be
			// allocated: compact BEFORE the first fill and the large
			// arrays are simply never created.
			sgrids.compact();
			sgrids.initialize_kinematics_with_ignore(irregular);
		}
		else {
			// The caller wants the WHOLE grid at t = 0 (a solver
			// initialising its full domain needs that), so fill dense and
			// hand the memory back on the first time advance instead --
			// see lsGridSpline::update().
			sgrids.initialize_kinematics(irregular);
		}
	}

	// ---- swd ------------------------------------------------------
	if (inputdata.wavetype >= 30 && inputdata.wavetype < 40) {
#if defined(SWD_enable)
		cfdwm_mpi::out() << inputdata.swdFileName.c_str() << std::endl;
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

		cfdwm_mpi::out() << "************************************************************************* " << std::endl;
		cfdwm_mpi::out() << "The following info is stored in the specified swd file " << inputdata.swdFileName << std::endl;
		cfdwm_mpi::out() << "************************************************************************* " << std::endl;
		for (size_t i = 0; i < v.size(); i++)
			cfdwm_mpi::out() << v[i] << std::endl;
		cfdwm_mpi::out() << "************************************************************************* " << std::endl;

		// Compute spectral wave properties from the SWD file (long-crested shapes only)
		inputdata.swd_spectral = compute_swd_spectral_properties(swd, inputdata.gravity, inputdata.depth);
		if (inputdata.swd_spectral.valid) {
			cfdwm_mpi::out() << "Spectral properties derived from SWD file:" << std::endl;
			cfdwm_mpi::out() << "  Hs = " << inputdata.swd_spectral.Hs << " m" << std::endl;
			cfdwm_mpi::out() << "  Tp = " << inputdata.swd_spectral.Tp << " s (peak period)" << std::endl;
			cfdwm_mpi::out() << "  T1 = " << inputdata.swd_spectral.T1 << " s (first-moment period, m0/m1)" << std::endl;
			cfdwm_mpi::out() << "  T2 = " << inputdata.swd_spectral.T2 << " s (zero-crossing period, sqrt(m0/m2))" << std::endl;
			cfdwm_mpi::out() << "************************************************************************* " << std::endl;
		}
#else
		cfdwm_mpi::err() << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
#endif
	}

	// SWD + lsgrid (case 34): sample the swd kinematics onto the spline grid.
	// Must come after the swd object is created above. The fill is serial
	// within the process (the swd library allows only one object per process
	// and is stateful); MPI ranks each fill their own column block.
	if (inputdata.wavetype == 34) {
#if defined(SWD_enable)
		if (!inputdata.lsgrid_domain_given && !inputdata.lsgrid_domain_override)
			fatal("INPUTFILE ERROR: the [lsgrid] extent is not defined. Give "
				"it either as\n    bounds <xmin> <xmax> <ymin> <ymax>\n"
				"in the [lsgrid] section, or by calling\n"
				"    wave_SetLsgridBounds(xmin, xmax, ymin, ymax)\n"
				"from the linked solver BEFORE wave_Initialize(). The latter "
				"is preferable when the solver already knows its own domain, "
				"since the two cannot then drift apart.");
		if (inputdata.lsgrid_domain_override)
			for (int b = 0; b < 4; b++)
				inputdata.lsgrid_domain[b] = inputdata.lsgrid_domain_set[b];
#if defined(SWD_enable)
		{
			// The SWD file states its own spectral resolution, so the grid
			// can be sized from it exactly as it is from our own component
			// set: the highest wavenumber carried is n*dk. Names follow the
			// SWD standard; a shape that does not define them throws, and we
			// fall back on the spectral estimate from the elevation DFT.
			double Lmin = 0., Tmin = 0.;
			bool have = false;
			try {
				const double dk = swd->GetReal("dk");
				const int    n  = swd->GetInt("n");
				if (dk > 0. && n > 0) {
					Lmin = 2.0 * M_PI / (n * dk);
					const double om = std::sqrt(inputdata.gravity * (n * dk));
					Tmin = 2.0 * M_PI / om;
					have = true;
				}
			}
			catch (...) { swd->ExceptionClear(); }
			if (!have && inputdata.swd_spectral.valid
					&& inputdata.swd_spectral.T2 > 0.) {
				// no usable metadata: fall back on the zero-crossing period
				Tmin = inputdata.swd_spectral.T2;
				const double om = 2.0 * M_PI / Tmin;
				Lmin = 2.0 * M_PI / (om * om / inputdata.gravity);
				have = true;
			}
			if (have)
				lsgrid_autosize(Lmin, Tmin, "SWD spectral resolution");
		}
#endif
		if (inputdata.lsgrid_border_override) {
			for (int b = 0; b < 4; b++)
				inputdata.lsgrid_border[b] = inputdata.lsgrid_border_set[b];
			inputdata.lsgrid_border_mode = true;
		}
		sgrids.domain[0] = inputdata.lsgrid_domain[0];
		sgrids.domain[1] = inputdata.lsgrid_domain[1];
		sgrids.domain[2] = inputdata.lsgrid_domain[2];
		sgrids.domain[3] = inputdata.lsgrid_domain[3];
		sgrids.nx = inputdata.lsgrid_nx;
		sgrids.ny = inputdata.lsgrid_ny;
		sgrids.nl = inputdata.lsgrid_nl;
		sgrids.t0 = inputdata.lsgrid_t0;
		sgrids.dt = inputdata.lsgrid_dt;
		sgrids.stretch_beta = inputdata.lsgrid_stretch_beta;
		// automatic stretching needs the discrete component set, which the SWD
		// backend does not expose; fall back to the fixed value
		sgrids.stretch_auto = false;
		if (inputdata.lsgrid_stretch_auto) {
			cfdwm_mpi::out() << "Note: 'stretch_beta auto' is not available with SWD "
				"input; using beta=" << sgrids.stretch_beta << std::endl;
		}
		// same for the automatic layer count
		sgrids.nl_auto = false;
		if (!inputdata.lsgrid_nl_set) {
			cfdwm_mpi::out() << "Note: automatic [lsgrid] nl is not available "
				"with SWD input; using nl=" << sgrids.nl << std::endl;
		}
		for (int b = 0; b < 4; b++) sgrids.border[b] = inputdata.lsgrid_border[b];
		sgrids.border_margin = inputdata.lsgrid_border_margin;
		sgrids.border_mode = inputdata.lsgrid_border_mode;
		sgrids.init_only = inputdata.lsgrid_init_only;
		sgrids.ignore_at_init = inputdata.lsgrid_ignore_at_init;
		// the SWD sampler does not fill a pressure field; wave_nhPres()
		// evaluates p_nh directly from the SWD object instead
		sgrids.store_pressure = false;
		if (inputdata.lsgrid_store_pressure == 1)
			cfdwm_mpi::out() << "Note: [lsgrid] 'store_pressure' is not "
				"available with SWD input; wave_nhPres() is served directly "
				"from the SWD field." << std::endl;
		sgrids.allocate();

		sgrids.water_depth = inputdata.depth;
		sgrids.swl = inputdata.swl;
		sgrids.set_ignore();

		if (sgrids.ignore_at_init) {
			sgrids.initialize_kinematics_with_ignore(swd);
		}
		else {
			sgrids.initialize_kinematics(swd);
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
	if (inputdata.wavetype == 51) {
		double yaw = inputdata.slope_direction * PI / 180.;
		double pitch = inputdata.slope_angle * PI / 180.;
		// Sloping water at rest. Compute rotation matrix from given angles for slope
		inputdata.slope_rotmat[0] = cos(yaw) * cos(pitch) * cos(yaw) + sin(yaw) * sin(yaw);
		inputdata.slope_rotmat[1] = cos(yaw) * cos(pitch) * sin(yaw) - sin(yaw) * cos(yaw);
		inputdata.slope_rotmat[2] = -cos(yaw) * sin(pitch);
		inputdata.slope_rotmat[3] = sin(yaw) * cos(pitch) * cos(yaw) - cos(yaw) * sin(yaw);
		inputdata.slope_rotmat[4] = sin(yaw) * cos(pitch) * sin(yaw) + cos(yaw) * cos(yaw);
		inputdata.slope_rotmat[5] = -sin(yaw) * sin(pitch);
		inputdata.slope_rotmat[6] = sin(pitch) * cos(yaw);
		inputdata.slope_rotmat[7] = sin(pitch) * sin(yaw);
		inputdata.slope_rotmat[8] = cos(pitch);
	}
	cfdwm_mpi::out() << "WaveID: " << inputdata.wavetype << std::endl;

	return 0;
}

} // namespace

// ------------------------------------------------------- entry point

int cfdwm_read_input() {
	// Re-initialization support: a second wave_Initialize() in the same
	// process must start from a clean parse state, or the section
	// uniqueness flags left by the previous run trigger a spurious
	// "duplicate section" error. Reset everything to defaults, carrying
	// across only the programmatic pre-init overrides
	// (wave_SetLsgridBounds/Borders/SetBorderMargin, which are
	// documented to be called before initialization).
	{
		CFDwavemakerInputdata fresh;
		for (int i = 0; i < 4; ++i) {
			fresh.lsgrid_domain_set[i] = inputdata.lsgrid_domain_set[i];
			fresh.lsgrid_border_set[i] = inputdata.lsgrid_border_set[i];
		}
		fresh.lsgrid_domain_override = inputdata.lsgrid_domain_override;
		fresh.lsgrid_border_override = inputdata.lsgrid_border_override;
		fresh.lsgrid_border_margin   = inputdata.lsgrid_border_margin;
		inputdata = fresh;
	}
	// The parse phase also writes directly into the wave-object
	// singletons (extrapolation method, ramps, component tables, ...),
	// so they must be re-constructed as well or settings leak from the
	// previous initialization into the next one. Placement new is used
	// because Irregular contains a mutex and is therefore not
	// assignable; the same pattern is applied to all of them for
	// uniformity.
	irregular.~Irregular();   new (&irregular) Irregular();
	stokes5.~Stokes5();       new (&stokes5) Stokes5();
	fenton.~FentonStream();   new (&fenton) FentonStream();
	wavemaker.~Wavemaker();   new (&wavemaker) Wavemaker();
	sgrids.~lsGridSpline();   new (&sgrids) lsGridSpline();
	ramp.~Ramp();             new (&ramp) Ramp();
	std::string lineA;
	std::ifstream fid;
	std::string res;

	// Check for waveinput.dat file in most common locations
	fid.open("./waveinput.dat");
	// check one step up in the folder tree (this is used in the latest comflow version)
	if (fid.fail()) {
		fid.open("./constant/waveinput.dat");
		if (!fid.fail()) {
			cfdwm_mpi::out() << "Waveinput.dat file found in folder ./constant/" << std::endl;
		}
	}
	// Special cases for comflow
	if (fid.fail()) {
		fid.open("../input_files/waveinput.dat");
		if (!fid.fail()) {
			cfdwm_mpi::out() << "Waveinput.dat file found in folder ./input_files/" << std::endl;
		}
	}
	if (fid.fail()) {
		fid.open("../waveinput.dat");
		if (!fid.fail()) {
			cfdwm_mpi::out() << "Waveinput.dat file found in folder ../" << std::endl;
		}
	}
	if (fid.fail()) {
		// printed on every rank: a missing file on one node only is exactly
		// the case that must not be hidden by master-only logging
		if (cfdwm_mpi::size() > 1)
			std::cerr << "[rank " << cfdwm_mpi::rank() << "] ";
		std::cerr << "Could not open waveinput.dat (is it really there?) " << std::endl;
		return -1;
	}
	// READ INPUT FILE AND REMOVE COMMENT LINES
	while (fid.good()) {
		getline(fid, lineA);
		lineA.erase(std::find(lineA.begin(), lineA.end(), '#'), lineA.end());
		if (lineA.length() > 0) {
			res += lineA + "\n";
		}
	}
	fid.close();

	std::istringstream f(res);

	// check file version (first line only)
	getline(f, lineA);
	trim(lineA);
	if (lineA.substr(0, 2).compare("@v")) {
		fatal("Please specify the input file version on the first line of "
			"the input file (hint: @v310)");
	}
	if (stoi(lineA.substr(2, 3)) < 300) {
		fatal(std::string("your waveinput.dat file is outdated and does not "
			"match the current version of CFDwavemaker (")
			+ CFDWM_VERSION + "). See the manual for the updated format.");
	}

	cfdwm_mpi::out() << std::boolalpha; // display booleans as true or false when printed (instead of 0 or 1)

	parse_all_sections(f);
	return resolve_and_initialize();
}
