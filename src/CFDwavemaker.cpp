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
#include "layer_stretch.h"
#include "sigma_fft.h"
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


// Variables
//int nfreq, ndir, wavetype, extmet, pertmet, meth, bandwidth, n_timesteps, rampswitch, normalizeA, spreadfunc;




#include "input_reader.h"

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

// The waveinput.dat reader (two-phase: parse -> resolve/initialize)
// lives in input_reader.cpp; see input_reader.h.

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
	// below-mudline guard: with variable bathymetry ([seabed]) the
	// check is against the LOCAL bed level h(x, y)
	if (irregular.seabed.active())
		zpt = std::max(inputdata.swl - irregular.seabed.depth(xpt, ypt), zpt);
	else
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
		cfdwm_mpi::err() << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}


	default:
		return temp;
	}
}


double wave_VeloY(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	// below-mudline guard: with variable bathymetry ([seabed]) the
	// check is against the LOCAL bed level h(x, y)
	if (irregular.seabed.active())
		zpt = std::max(inputdata.swl - irregular.seabed.depth(xpt, ypt), zpt);
	else
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
		cfdwm_mpi::err() << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		return 0.0;
	}
}


double wave_VeloZ(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	// below-mudline guard: with variable bathymetry ([seabed]) the
	// check is against the LOCAL bed level h(x, y)
	if (irregular.seabed.active())
		zpt = std::max(inputdata.swl - irregular.seabed.depth(xpt, ypt), zpt);
	else
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
		cfdwm_mpi::err() << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		return 0.0;
	}
}

// todo: not fully functional for all wavetypes, but not really used either....
double wave_DynPres(double xpt, double ypt, double zpt, double tpt)
{
	// Quickfix 07022018 To avoid issues with values below mudline
	// below-mudline guard: with variable bathymetry ([seabed]) the
	// check is against the LOCAL bed level h(x, y)
	if (irregular.seabed.active())
		zpt = std::max(inputdata.swl - irregular.seabed.depth(xpt, ypt), zpt);
	else
		zpt = std::max(-inputdata.depth, zpt);

	switch (inputdata.wavetype) {
		// irregular waves
	case 1:
		return ramp.ramp(tpt, xpt, ypt) * irregular.dp(tpt, xpt, ypt, zpt);
	// irregular + lsgrid: interpolated from the stored pressure field when
	// available (v3.7.12), else evaluated pointwise (slow but correct)
	case 5:
	{
		if (sgrids.store_pressure) {
			std::vector<double> v = sgrids.get_kinematics_at_point(
				tpt, xpt, ypt, zpt, inputdata.depth, true);
			return ramp.ramp(tpt, xpt, ypt) * v[4];
		}
		return ramp.ramp(tpt, xpt, ypt) * irregular.dp(tpt, xpt, ypt, zpt);
	}

	default:
		return 0.0;
	}
}

// Nonhydrostatic pressure: p_nh = p_total - rho*g*(eta - z)
double wave_nhPres(double xpt, double ypt, double zpt, double tpt)
{
	// Clamp to seabed
	// below-mudline guard: with variable bathymetry ([seabed]) the
	// check is against the LOCAL bed level h(x, y)
	if (irregular.seabed.active())
		zpt = std::max(inputdata.swl - irregular.seabed.depth(xpt, ypt), zpt);
	else
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
	// irregular + lsgrid (wavetype 5). With store_pressure (the default
	// for border-mode grids since v3.7.12) the grid carries the full
	// second-order dynamic pressure as a spline field -- filled by the
	// sigma-FFT fill or the pointwise fill alongside the velocities --
	// and p_nh is interpolated at grid speed with the same truncation as
	// the served kinematics.
	// Without it, fall back to the FIRST-ORDER p_nh by direct summation
	// (v3.7.11 behaviour): O(nfreq) per call, consistent truncation
	// (dp1 - rho*g*eta1). Before v3.7.11 this silently returned 0 --
	// generating-absorbing boundary conditions that subtract the incident
	// non-hydrostatic pressure then over-generate (measured ~+20%
	// amplitude on a directional coastal sea state).
	case 5:
	{
		if (sgrids.store_pressure) {
			std::vector<double> v = sgrids.get_kinematics_at_point(
				tpt, xpt, ypt, zpt, inputdata.depth, true);
			// v[0] is the interpolated elevation INCLUDING swl -- the same
			// convention as eta_raw in case 1 above
			return ramp.ramp(tpt, xpt, ypt)
				* (v[4] - irregular.rho * irregular.g * v[0]);
		}
		double eta1_raw = irregular.eta1(tpt, xpt, ypt);
		double dp1_raw = (irregular.vertical_theory == 1)
			? irregular.dp1_smit(tpt, xpt, ypt, zpt)
			: irregular.dp1(tpt, xpt, ypt, std::min(0., zpt - inputdata.swl));
		return ramp.ramp(tpt, xpt, ypt) * (dp1_raw - irregular.rho * irregular.g * eta1_raw);
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
		return 0.0;
#endif
	}
	// SWD - computed from primitives: p_nh = -rho*(dphi/dt + 0.5*|u|^2 + g*eta)
	// (case 34, swd + lsgrid: the spline grid stores no pressure field, so
	// p_nh is evaluated directly from the SWD object -- stateful-serial but
	// exact. Before v3.7.11 case 34 silently returned 0 here.)
	case 34:
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
		cfdwm_mpi::err() << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of swd library specified in the waveinput.dat file. Please recompile CFDwavemaker with SWD_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
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
		cfdwm_mpi::err() << "Use of vtk library specified in the waveinput.dat file. Please recompile CFDwavemaker with VTK_enable=1." << std::endl;
#endif
	}
	default:
		// variable bathymetry: local bed level from the [seabed] surface
		// (min-depth clamped, i.e. what the kinematics would see)
		if (irregular.seabed.loaded())
			return inputdata.swl - irregular.seabed.depth(xpt, ypt);
		return -inputdata.depth;;
	}
}

double wave_SeabedRaw(double xpt, double ypt)
{
	// True (unclamped) bed level from the [seabed] surface: no
	// minimum-depth limit, so land is returned as positive elevation.
	// For solvers importing the real bathymetry (wetting/drying, run-up);
	// the wave kinematics themselves always use the clamped depth.
	if (irregular.seabed.loaded())
		return inputdata.swl - irregular.seabed.depth_raw(xpt, ypt);
	return -inputdata.depth;
}

void wave_SeabedRelease(void)
{
	// Free the [seabed] mesh once init-time sampling is done: in
	// geometry-only mode (kinematics 0) the surface is only queried
	// while a linked solver imports its bathymetry (zb, phase-speed and
	// layer-distribution precomputes), yet a large OBJ costs hundreds
	// of MB per rank for the whole run. Refused (with a warning) when
	// the kinematics use the surface; queries after release abort.
	irregular.seabed.release();
}

int wave_SeabedActive(void)
{
	// Lets a linked solver distinguish a real [seabed] surface from the
	// constant-depth fallback of wave_Seabed() before sampling bathymetry.
	// True also in geometry-only mode ([seabed] kinematics 0).
	return irregular.seabed.loaded() ? 1 : 0;
}

int wave_SetLsgridBounds(double xmin, double xmax, double ymin, double ymax)
{
	// A linked solver knows its own domain; let it say so rather than
	// duplicating the numbers in waveinput.dat where they can drift apart.
	// Must precede wave_Initialize(): the grid is built during init.
	if (CFDwavemaker_is_initialized()) {
		cfdwm_mpi::err() << "wave_SetLsgridBounds(): ignored, must be called "
			"BEFORE wave_Initialize()." << std::endl;
		return 0;
	}
	inputdata.lsgrid_domain_set[0] = xmin;
	inputdata.lsgrid_domain_set[1] = xmax;
	inputdata.lsgrid_domain_set[2] = ymin;
	inputdata.lsgrid_domain_set[3] = ymax;
	inputdata.lsgrid_domain_override = true;
	return 1;
}

int wave_SetBorderMargin(double metres)
{
	// How far inward the client will query. Without it the stored strip is
	// sized in NODES, which silently stops being enough when the client's
	// cells are coarse relative to the grid spacing (a coarse CFD grid with
	// a finely resolved wave lattice). Given it, the strip is sized in
	// metres and the mismatch cannot occur. Must precede wave_Initialize().
	if (CFDwavemaker_is_initialized()) {
		cfdwm_mpi::err() << "wave_SetBorderMargin(): ignored, must be called "
			"BEFORE wave_Initialize()." << std::endl;
		return 0;
	}
	inputdata.lsgrid_border_margin = metres > 0. ? metres : 0.;
	return 1;
}

int wave_SetLsgridBorders(int left, int right, int bottom, int top)
{
	// Which borders actually need kinematics. Everything else is dropped
	// from storage. This is the resolution-independent way to say it --
	// unlike ignore_subdomain, which describes the complement in absolute
	// metres and has to be re-derived whenever the client's grid changes.
	// Must precede wave_Initialize().
	if (CFDwavemaker_is_initialized()) {
		cfdwm_mpi::err() << "wave_SetLsgridBorders(): ignored, must be called "
			"BEFORE wave_Initialize()." << std::endl;
		return 0;
	}
	inputdata.lsgrid_border_set[0] = left  ? 1 : 0;
	inputdata.lsgrid_border_set[1] = right ? 1 : 0;
	inputdata.lsgrid_border_set[2] = bottom ? 1 : 0;
	inputdata.lsgrid_border_set[3] = top   ? 1 : 0;
	inputdata.lsgrid_border_override = true;
	return 1;
}

int wave_SeabedBounds(double* xmin, double* xmax, double* ymin, double* ymax)
{
	// xy bounding box of the [seabed] footprint, so a linked solver can
	// check its domain against the bathymetry BEFORE it starts. Returns 0
	// (and leaves the outputs untouched) when no surface is loaded.
	// NB: passing this test is necessary, not sufficient -- the mesh need
	// not fill its own bbox. Use wave_SeabedContains() for the exact test.
	if (!irregular.seabed.loaded())
		return 0;
	double a, b, c, d;
	irregular.seabed.bounds(a, b, c, d);
	if (xmin) *xmin = a;
	if (xmax) *xmax = b;
	if (ymin) *ymin = c;
	if (ymax) *ymax = d;
	return 1;
}

int wave_SeabedContains(double xpt, double ypt)
{
	// Exact point-in-footprint test: is (x, y) covered by a triangle of
	// the [seabed] surface? No clamping and no warning, unlike a depth
	// query. 0 when outside, or when no surface is loaded.
	return irregular.seabed.contains(xpt, ypt) ? 1 : 0;
}

long wave_SeabedOutsideCount(void)
{
	// How many depth queries have been served by out-of-footprint
	// clamping so far. Non-zero means part of the sampled region is not
	// actually described by the bathymetry: the bed there is the flat,
	// zero-slope skirt of the nearest mesh point, which is a common cause
	// of boundary blow-up in a linked solver.
	return irregular.seabed.loaded() ? irregular.seabed.outside_count() : 0;
}

double wave_LocalDepth(double xpt, double ypt)
{
	// Local still-water depth [m], positive down. Whenever a [seabed]
	// surface is loaded it is THE depth, whatever generates the kinematics
	// (direct summation, lsgrid / sigma-FFT, ...) and whether or not the
	// surface feeds the kinematics themselves -- geometry-only mode
	// ([seabed] kinematics 0) still reports the real bathymetry here.
	// Without a [seabed] surface the constant reference depth from the
	// [general input data] 'depth' line is returned instead.
	// Unclamped by min_depth, so land comes back negative: consistent with
	// wave_SeabedRaw(), which is what a solver imports as its bathymetry.
	if (irregular.seabed.loaded())
		return irregular.seabed.depth_raw(xpt, ypt);
	return inputdata.depth;
}

/* ---------------------------------------------------------------------------
   Optimal geometric layer grading for a linked layered solver.

   The grading minimises the vertical reconstruction error on the sea state's
   spectral RMS velocity profile (layer_stretch.h, the same optimiser the
   lsgrid uses for its sigma levels) and is returned as the TOP-LAYER-TO-
   UNIFORM thickness ratio in (0, 1] -- the parameter consumed directly by
   e.g. Basilisk's geometric_beta() (multilayer 'cell_aspect_ratio').
   1.0 = uniform; negative on failure.

   optimal_layer_ratio_at() evaluates it at an ARBITRARY depth so a solver
   with variable bathymetry can give every column its own layer distribution
   rather than stretching one deep-water distribution over the whole domain.
   The component set is carried to the local depth by conserving each
   component's frequency -- its wavenumber then follows the local dispersion
   relation -- and its energy flux (Green's law for the amplitude). At the
   reference depth this reproduces the stored components exactly. Refraction
   and directional redistribution are ignored: only the SHAPE of the vertical
   velocity profile matters here, and they barely move it.
   --------------------------------------------------------------------------- */

// Newton on g k tanh(k h) = om^2, starting from the reference-depth k.
static double local_wavenumber(double om, double h, double g, double k0)
{
	double k = (k0 > 0.) ? k0 : om * om / g;
	for (int it = 0; it < 50; ++it) {
		const double th = std::tanh(k * h);
		const double dk = (g * k * th - om * om)
			/ (g * (th + k * h * (1. - th * th)));
		k -= dk;
		if (!(k > 0.))				// Newton overshoot in very shallow water
			k = om / std::sqrt(g * h);	// shallow-water restart
		if (std::fabs(dk) < 1e-12 * k)
			break;
	}
	return k;
}

// cg = (om/2k)(1 + 2kh/sinh 2kh), written to survive both limits of kh
static double group_velocity(double om, double k, double h)
{
	if (!(k > 0.))
		return 0.;
	const double kh2 = 2. * k * h;
	const double f = (kh2 > 700.) ? 0. : (kh2 < 1e-6 ? 1. : kh2 / std::sinh(kh2));
	return 0.5 * om * (1. + f) / k;
}

static double optimal_layer_ratio_at(int nlayers, double depth)
{
	const int N = irregular.nfreq * irregular.ndir;
	if (nlayers < 3 || N <= 0 || (int)irregular.k.size() < N
			|| (int)irregular.omega.size() < N
			|| (int)irregular.A.size() < N
			|| depth <= 0. || inputdata.depth <= 0.)
		return -1.;
	const double g = inputdata.gravity, href = inputdata.depth;
	std::vector<double> c(N), k(N);
	for (int j = 0; j < N; ++j) {
		const double om = irregular.omega[j];
		if (!(om > 0.)) {
			c[j] = 0.;
			k[j] = irregular.k[j];
			continue;
		}
		const double kloc = local_wavenumber(om, depth, g, irregular.k[j]);
		const double cg_ref = group_velocity(om, irregular.k[j], href);
		const double cg_loc = group_velocity(om, kloc, depth);
		// Green's law: A ~ cg^{-1/2}
		const double sh = (cg_loc > 0. && cg_ref > 0.)
			? std::sqrt(cg_ref / cg_loc) : 1.;
		// surface horizontal velocity amplitude at the local depth
		const double v = irregular.A[j] * sh * g * kloc / om;
		c[j] = v * v;
		k[j] = kloc;
	}
	// nlayers layers are bounded by nlayers+1 sigma levels
	const double beta = layer_stretch::optimal_beta(nlayers + 1, c, k, depth);
	if (beta < 0.)
		return -1.;
	// adjacent-layer thickness ratio of the exponential mapping, then
	// inverted through geometric_beta's r(rmin) parametrisation
	const double r = std::exp(beta / double(nlayers));
	return 1.0 / (1.0 + 0.5 * (r - 1.0) * double(nlayers - 1));
}

/* One evaluation costs ~10^7 flops (the optimiser sweeps ~180 candidate beta,
   each rebuilding the profile over the whole component set), far too much to
   call per grid column. The ratio is smooth and monotone in depth, so it is
   tabulated once -- log-spaced, since the interesting variation is in the
   first few metres -- and interpolated afterwards. */

namespace {
	const double RATIO_TABLE_HMIN = 0.05;	// [m], shallow end of the table
	const int RATIO_TABLE_N = 48;

	struct RatioTable {
		bool built = false;	// build attempted for `nlayers` (success or not)
		bool valid = false;
		int nlayers = -1;
		double lnh0 = 0., lnh1 = 0.;	// uniform in ln(depth): no search needed
		std::vector<double> ratio;
	};
	RatioTable ratio_table;
}

static void build_ratio_table(int nlayers)
{
	RatioTable t;
	t.built = true;
	t.nlayers = nlayers;
	double hmax = inputdata.depth;
	if (irregular.seabed.loaded())
		hmax = std::max(hmax, irregular.seabed.max_depth());
	if (!(hmax > RATIO_TABLE_HMIN))
		hmax = 10. * RATIO_TABLE_HMIN;
	hmax *= 1.05;			// headroom for surge / a deeper solver bed
	t.lnh0 = std::log(RATIO_TABLE_HMIN);
	t.lnh1 = std::log(hmax);
	t.ratio.resize(RATIO_TABLE_N);
	for (int i = 0; i < RATIO_TABLE_N; ++i) {
		const double lh = t.lnh0
			+ (t.lnh1 - t.lnh0) * double(i) / double(RATIO_TABLE_N - 1);
		const double r = optimal_layer_ratio_at(nlayers, std::exp(lh));
		// a failure at one depth is not a reason to lose the whole table:
		// the shallow end degenerates to uniform layers anyway
		t.ratio[i] = (r > 0.) ? r : 1.0;
		if (r > 0.)
			t.valid = true;
	}
	ratio_table = t;
}

double wave_OptimalLayerRatioAtDepth(int nlayers, double depth)
{
	if (nlayers < 3)
		return -1.;
	// The table is built once and then read-only. Callers that evaluate it
	// from a parallel region should prime it with one call beforehand; the
	// critical section below is the belt-and-braces for those that do not.
	if (!(ratio_table.built && ratio_table.nlayers == nlayers)) {
#if defined(_OPENMP)
#pragma omp critical(cfdwm_ratio_table)
#endif
		{
			if (!(ratio_table.built && ratio_table.nlayers == nlayers))
				build_ratio_table(nlayers);
		}
	}
	if (!ratio_table.valid)
		return -1.;
	if (!(depth > 0.))
		return 1.0;		// dry / land: uniform layers
	const double lh = std::log(depth);
	if (lh <= ratio_table.lnh0)
		return ratio_table.ratio.front();
	if (lh >= ratio_table.lnh1)
		return ratio_table.ratio.back();
	const double s = (lh - ratio_table.lnh0)
		/ (ratio_table.lnh1 - ratio_table.lnh0) * double(RATIO_TABLE_N - 1);
	int i = int(s);
	if (i > RATIO_TABLE_N - 2)
		i = RATIO_TABLE_N - 2;
	const double w = s - double(i);
	return ratio_table.ratio[i] * (1. - w) + ratio_table.ratio[i + 1] * w;
}

double wave_OptimalLayerRatioAt(int nlayers, double xpt, double ypt)
{
	return wave_OptimalLayerRatioAtDepth(nlayers, wave_LocalDepth(xpt, ypt));
}

double wave_OptimalLayerRatio(int nlayers)
{
	// the reference ('depth' line) depth, evaluated exactly
	return optimal_layer_ratio_at(nlayers, inputdata.depth);
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
	// MPI: every rank sees the same wave field; only rank 0 appends to
	// the probe file (concurrent appends from all ranks would garble it)
	if (cfdwm_mpi::rank() != 0) return;
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
	// MPI: only the master rank echoes the banner and the input summary
	// (every rank reads the same waveinput.dat, and N copies of the full
	// summary makes the job log unusable). See cfdwm_mpi::out().
	cfdwm_mpi::out() << "\n\n***********************************************\n\n" << std::endl;
	cfdwm_mpi::out() << "---------------------------------------" << std::endl;
	cfdwm_mpi::out() << "CFD WAVEMAKER " << CFDWM_VERSION << std::endl;
	cfdwm_mpi::out() << "---------------------------------------" << std::endl;

	std::string cwd = get_current_dir();
	cfdwm_mpi::out() << "working directory: " << cwd << std::endl;
	if (cfdwm_mpi::size() > 1)
		cfdwm_mpi::out() << "MPI: " << cfdwm_mpi::size()
			<< " ranks (log written by rank 0 only)" << std::endl;

	if (cfdwm_read_input() != 0) {
		return -1;
	}
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
	// drop the layer-grading table: a re-initialize may bring a different
	// spectrum or bathymetry
	ratio_table = RatioTable();
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
	case 5: // lsgrid: same irregular component set
	case 1:
		return irregular.phase_speed(opt);
	default:
		return 0.;
	}
}

double wave_PeakPeriod(void)
{
	// spectral peak period Tp: the largest-amplitude component (same
	// estimate as the init summary print)
	if ((inputdata.wavetype == 1 || inputdata.wavetype == 5)
			&& !irregular.A.empty()) {
		int imax = 0;
		for (int i = 1; i < (int)irregular.A.size(); ++i)
			if (irregular.A[i] > irregular.A[imax])
				imax = i;
		return 2.*M_PI/irregular.omega[imax];
	}
	return 0.;
}

double wave_PhaseVelocity(double xpt, double ypt, double T0)
{
	// Linear phase speed c = omega0/k at the LOCAL water depth, for a
	// fixed representative period T0 (T0 <= 0: spectral mean period T1,
	// falling back to the peak period -- the mean period is the more
	// representative "typical wave" scale and is slightly shorter).
	// The depth is wave_LocalDepth(): the [seabed] surface when one is
	// loaded (also in geometry-only mode, unclamped -- consistent with what
	// a linked solver imports via wave_SeabedRaw), else the constant
	// reference depth. Returns 0 on dry/land points or without a valid
	// period. Intended for once-per-cell precomputation of shoaling-
	// aware normalisations (e.g. kinematic breaking indicators u_s/c).
	if (T0 <= 0.)
		T0 = wave_mean_period(1);
	if (T0 <= 0.)
		T0 = wave_PeakPeriod();
	if (T0 <= 0.)
		return 0.;
	const double h = wave_LocalDepth(xpt, ypt);
	if (h <= 0.)
		return 0.;
	const double om = 2.*M_PI/T0;
	return om/local_wavenumber(om, h, inputdata.gravity, 0.);
}

double wave_mean_length(int opt) {
	switch (inputdata.wavetype) {
		// Linear wave theory, expenential profile used above free surface
	case 5: // lsgrid: same irregular component set
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
	case 5: // lsgrid: same irregular component set
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