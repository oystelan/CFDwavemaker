#ifndef Utils_H
#define Utils_H

#include "Irregular.h"
#include <sstream>

#if defined(SWD_enable)
#include "SpectralWaveData.h"

// Spectral wave properties derived from an SWD object by sampling the surface
// elevation along one periodic length in space and recovering the spectral
// amplitudes via a DFT. Currently supports shape 1 (long-crested) only.
struct SwdSpectralProps {
    double m0 = 0.0;   // 0th spectral moment (variance)
    double m1 = 0.0;   // 1st spectral moment
    double m2 = 0.0;   // 2nd spectral moment
    double Hs = 0.0;   // Significant wave height = 4*sqrt(m0)
    double Tp = 0.0;   // Peak period
    double T1 = 0.0;   // First-moment period = 2*pi*m0/m1
    double T2 = 0.0;   // Zero-crossing period = 2*pi*sqrt(m0/m2)
    bool valid = false;
};

// Compute spectral properties from an initialized SWD object.
// `gravity` and `depth` are used for the linear dispersion relation
// omega^2 = g*k*tanh(k*h). Pass depth <= 0 for deep water.
SwdSpectralProps compute_swd_spectral_properties(SpectralWaveData* swd,
                                                 double gravity, double depth);
#endif

// Ramp class contains various ramp functions
class Ramp {
private:
	double ramp1d(double x, double xstart, double xend, bool inv);
public:
	int timeramp, xramp, yramp; // logical operators for switching on and of rampfunctions
	double time_rampup_start, time_rampup_end, time_rampdown_start, time_rampdown_end;
	double x_rampup_start, x_rampup_end, x_rampdown_start, x_rampdown_end;
	double y_rampup_start, y_rampup_end, y_rampdown_start, y_rampdown_end;
	bool ramp_init = false;
	bool ramp_init_time_up = false;
	bool ramp_init_time_down = false;
	bool ramp_init_x_up = false;
	bool ramp_init_x_down = false;
	bool ramp_init_y_up = false;
	bool ramp_init_y_down = false;

	double ramp(double time, double x, double y);
	
};

int dirExists(const char* const path);

void createDirectory(std::string sPath);
std::string get_current_dir();

#endif
