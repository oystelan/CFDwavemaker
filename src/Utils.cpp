
#include "Utils.h"
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <sys/stat.h>
//#if defined(_WIN32)
//#include < direct.h >
//#endif
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>

// -------------------------------------------------------------------------------------------------
// ramp class function
// -------------------------------------------------------------------------------------------------

//Define some useful functions
/* Rampfunction */
// NB: Not yet implemented inverse ramp
double Ramp::ramp1d(double x, double xstart, double xend, bool inv) {

	if (inv){
		return std::max(std::min(1. - ((x - xstart) / (xend - xstart)),1.), 0.);
	}
	else {
		return std::max(std::min(((x - xstart) / (xend - xstart)), 1.), 0.);
	}
}

bool comp(double a, double b)
{
	return (a < b);
}

// Todo: implement ramps
double Ramp::ramp(double t, double x, double y) {
	if (ramp_init) {
		double ramps[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
		if (ramp_init_time_up) {
			ramps[0] = ramp1d(t, time_rampup_start, time_rampup_end, false);
		}
		if (ramp_init_time_down) {
			ramps[1] = ramp1d(t, time_rampdown_start, time_rampdown_end, true);
		}
		if (ramp_init_x_up) {
			ramps[2] = ramp1d(x, x_rampup_start, x_rampup_end, false);
		}
		if (ramp_init_x_down) {
			ramps[3] = ramp1d(x, x_rampdown_start, x_rampdown_end, true);
		}
		if (ramp_init_y_up) {
			ramps[4] = ramp1d(y, y_rampup_start, y_rampup_end, false);
		}
		if (ramp_init_y_down) {
			ramps[5] = ramp1d(y, y_rampdown_start, y_rampdown_end, true);
		}
		return *std::min_element(std::begin(ramps), std::end(ramps));
	}
	else {
		return 1.0;
	}
}


/******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full path AND if path leaf is a dir.
 *
 * @return  >0 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          <0 if an error occurred (errno is also set)
 *****************************************************************************/
int dirExists(const char* const path)
{
	struct stat info;

	int statRC = stat(path, &info);
	if (statRC != 0)
	{
		if (errno == ENOENT) { return 0; } // something along the path does not exist
		if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
		return -1;
	}

	return (info.st_mode & S_IFDIR) ? 1 : 0;
}

void createDirectory(std::string sPath) {


	int nError = 0;
#if defined(_WIN32)

	nError = _mkdir(sPath.c_str()); // can be used on Windows
#else
	mode_t nMode = 0733; // UNIX style permissions
	nError = mkdir(sPath.c_str(), nMode); // can be used on non-Windows
#endif
	if (nError != 0) {
		// handle your error here
	}
}


using namespace std;

std::string get_current_dir() {
	char buff[FILENAME_MAX]; //create string buffer to hold path
	GetCurrentDir(buff, FILENAME_MAX);
	string current_working_dir(buff);
	return current_working_dir;
}


#if defined(SWD_enable)

/*
Compute spectral wave properties from an SWD object.

Strategy: the spectral wavenumbers are deterministic (k_j = j*dk for shape 1),
so we recover the spectral amplitudes |h_j| by sampling eta(x) at N equispaced
points along one full periodic length L = 2*pi/dk and performing an explicit
DFT. For N >= 2*n_swd (Nyquist) the recovery is exact.

For shape 1, the elevation is:
    eta(x,t) = sum_j Re{ h_j(t) * exp(-i*k_j*x) }
             = sum_j [ Re(h_j) cos(k_j*x) + Im(h_j) sin(k_j*x) ]

With samples x_n = n*L/N, the inverse-DFT relations give the spectral
amplitudes A_j = |h_j| = sqrt(a_j^2 + b_j^2), with energy E_j = (1/2)*A_j^2.

Spectral moments: m_p = sum_j omega_j^p * E_j.
Mean periods: T1 = 2*pi*m0/m1, T2 = 2*pi*sqrt(m0/m2), Tp = 2*pi / omega_argmax.
*/
SwdSpectralProps compute_swd_spectral_properties(SpectralWaveData* swd,
                                                 double gravity, double depth)
{
    SwdSpectralProps p;

    if (swd == nullptr) {
        std::cerr << "compute_swd_spectral_properties: null SWD pointer." << std::endl;
        return p;
    }

    int shp = 0;
    try {
        shp = swd->GetInt("shp");
    } catch (...) {
        std::cerr << "compute_swd_spectral_properties: failed to query SWD 'shp'." << std::endl;
        return p;
    }

    if (shp != 1 && shp != 2) {
        std::cerr << "compute_swd_spectral_properties: only SWD shape 1 (deep) and 2 (finite depth) "
                  << "long-crested are currently supported. shp=" << shp << std::endl;
        return p;
    }

    // Use SWD's own depth if shape provides one (shape 2). Shape 1 reports d=-1 (deep water).
    double d = depth;
    try {
        double d_swd = swd->GetReal("d");
        if (d_swd > 0.0) d = d_swd;
    } catch (...) {
        // d remains the caller-supplied value
    }

    double dk = 0.0;
    int n_swd = 0, nsumx = -1;
    try {
        dk = swd->GetReal("dk");
        n_swd = swd->GetInt("n");
        nsumx = swd->GetInt("nsumx");
    } catch (...) {
        std::cerr << "compute_swd_spectral_properties: failed to query SWD parameters." << std::endl;
        return p;
    }

    if (nsumx > 0 && nsumx < n_swd) n_swd = nsumx;

    if (dk <= 0.0 || n_swd <= 0) {
        std::cerr << "compute_swd_spectral_properties: invalid dk or n (dk=" << dk
                  << ", n=" << n_swd << ")." << std::endl;
        return p;
    }

    const double L = 2.0 * M_PI / dk;
    const int N = 2 * n_swd;  // Nyquist sampling

    try {
        swd->UpdateTime(0.0);
    } catch (...) {
        std::cerr << "compute_swd_spectral_properties: UpdateTime(0) failed." << std::endl;
        return p;
    }

    // Sample eta at N equispaced x-positions over one period
    std::vector<double> eta_samples(N);
    for (int n = 0; n < N; ++n) {
        double xi = n * L / N;
        eta_samples[n] = swd->Elev(xi, 0.0);
    }

    // Explicit DFT for j = 1..n_swd (skip j=0 DC component)
    double tp_omega = 0.0;
    double max_amp2 = 0.0;

    for (int j = 1; j <= n_swd; ++j) {
        double a = 0.0, b = 0.0;
        for (int n = 0; n < N; ++n) {
            double phase = 2.0 * M_PI * j * n / N;
            a += eta_samples[n] * std::cos(phase);
            b += eta_samples[n] * std::sin(phase);
        }
        a *= 2.0 / N;
        b *= 2.0 / N;
        double amp2 = a * a + b * b;   // A_j^2

        double k_j = j * dk;
        double omega_j;
        if (d > 0.0) {
            omega_j = std::sqrt(gravity * k_j * std::tanh(k_j * d));
        } else {
            omega_j = std::sqrt(gravity * k_j);   // deep water
        }

        double E_j = 0.5 * amp2;
        p.m0 += E_j;
        p.m1 += omega_j * E_j;
        p.m2 += omega_j * omega_j * E_j;

        if (amp2 > max_amp2) {
            max_amp2 = amp2;
            tp_omega = omega_j;
        }
    }

    p.Hs = 4.0 * std::sqrt(std::max(p.m0, 0.0));
    p.Tp = (tp_omega > 0.0) ? 2.0 * M_PI / tp_omega : 0.0;
    p.T1 = (p.m1 > 0.0) ? 2.0 * M_PI * p.m0 / p.m1 : 0.0;
    p.T2 = (p.m2 > 0.0) ? 2.0 * M_PI * std::sqrt(p.m0 / p.m2) : 0.0;
    p.valid = true;
    return p;
}

#endif
